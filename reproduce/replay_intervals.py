"""Regenerate the interval benchmark and verify the two published result tables."""

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from functools import partial
import json
from math import isclose
from pathlib import Path

import numpy as np
from scipy.stats import beta, t

from crosslabel._field_kernel import (
    OUTER, block_se, coarse_sum, fine_sufficient, nested_estimate,
    spatial_estimate, value_gradient,
)
from interval_fields import field, merge_same_code, schedule


REFERENCE = Path(__file__).resolve().parent / "data" / "interval_qualification"
METHODS = {
    "nested": ("spatial_bootstrap_t_ci",),
    "primary": ("spatial_t_ci", "normal_ci"),
    "s3": ("studentized_ci", "percentile_ci", "delta_t_ci"),
    "s2": ("studentized_ci", "percentile_ci", "delta_t_ci"),
}


def coarse_intervals(x, weights, multiplicities):
    counts = x[:, :len(weights)].sum(0)
    if np.einsum("abc,a,b,c->", weights, counts, counts, counts) <= 0:
        return {}
    try:
        point, gradient = value_gradient(x.sum(0), weights)
        se = float(block_se(x, gradient))
        draws, gradients = value_gradient(multiplicities @ x, weights)
        se_star = block_se(x, gradients, multiplicities)
        low, high = np.percentile((draws - point) / se_star, [2.5, 97.5])
        critical = t.ppf(.975, len(x) - 1)
        return {
            "studentized_ci": [float(point - high * se), float(point - low * se)],
            "percentile_ci": np.percentile(draws, [2.5, 97.5]).tolist(),
            "delta_t_ci": [float(point - critical * se), float(point + critical * se)],
        }
    except (ValueError, FloatingPointError):
        return {}


def evaluate(xy, codes, radius, code_count, seed):
    fine, mapping, geometry = fine_sufficient(xy, codes, radius, code_count)
    blocks = geometry["coarse_blocks"]
    picks = np.random.default_rng(200000 + seed).integers(0, blocks, (OUTER, blocks))
    multiplicities = np.array([np.bincount(row, minlength=blocks) for row in picks])
    results = {}
    for mask, stats in fine.items():
        weights = stats["W"]
        coarse = coarse_sum(stats["x"], mapping, blocks)
        unadjusted = coarse_sum(stats["x_s2"], mapping, blocks)
        spatial = spatial_estimate(stats["x"], weights, geometry["shape"], blocks)
        nested, _ = nested_estimate(coarse, weights, spatial, seed)
        methods = dict(
            nested=nested, primary=spatial,
            s3=coarse_intervals(coarse, weights, multiplicities),
            s2=coarse_intervals(unadjusted, weights, multiplicities),
        )
        results[mask] = {
            mode: {name: methods[mode].get(name) for name in names}
            for mode, names in METHODS.items()
        }
    return results


def regenerate(plan):
    seed, condition, merge = plan
    xy, codes, truth, radius = field(seed, condition)
    views = [("native", xy, codes)]
    if merge:
        views.append(("merged_10nm", *merge_same_code(xy, codes)))
    code_count = 6 if condition.startswith("dense") else 4
    return [
        dict(seed=seed, condition=condition, representation=representation,
             truth=truth, records=len(points),
             intervals=evaluate(points, labels, radius, code_count, seed))
        for representation, points, labels in views
    ]


def rate_rule(usable, covered, positive, condition, planned):
    if usable != planned:
        return False
    low = 0. if covered == 0 else beta.ppf(.025, covered, usable - covered + 1)
    high = 1. if covered == usable else beta.ppf(.975, covered + 1, usable - covered)
    positives_ok = positive <= .05 * planned if condition.endswith("dimer") else positive >= .8 * planned
    return bool(covered >= np.ceil(.925 * planned) and low <= .95 <= high and positives_ok)


def covers(interval, truth):
    return interval[0] <= truth <= interval[1]


def outcomes(left, right):
    return {
        "both": sum(a and b for a, b in zip(left, right)),
        "nested_only": sum(a and not b for a, b in zip(left, right)),
        "comparator_only": sum(not a and b for a, b in zip(left, right)),
        "neither": sum(not a and not b for a, b in zip(left, right)),
    }


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["condition"], row["representation"]].append(row)
    summary, paired = [], []
    for (condition, representation), group in groups.items():
        planned = 400 if representation == "native" else 10
        for mask in ("all", "non_consecutive"):
            truths = [r["truth"] for r in group]
            methods = [r["intervals"][mask] for r in group]
            nested = [r["nested"]["spatial_bootstrap_t_ci"] for r in methods]
            for mode, names in METHODS.items():
                for name in names:
                    intervals = [r[mode][name] for r in methods]
                    valid = [(ci, y) for ci, y in zip(intervals, truths) if ci is not None]
                    covered = sum(covers(ci, y) for ci, y in valid)
                    positive = sum(ci[0] > 0 for ci, _ in valid)
                    identity = dict(condition=condition, representation=representation,
                                    mask=mask, mode=mode, method=name)
                    summary.append(dict(
                        **identity, planned=planned, observed=len(group), usable=len(valid),
                        unsupported=len(group) - len(valid), covered=covered, positive=positive,
                        mean_width=sum(ci[1] - ci[0] for ci, _ in valid) / len(valid) if valid else None,
                        original_rate_rule=rate_rule(len(valid), covered, positive, condition, planned)
                        if representation == "native" else None,
                    ))
                    if mode != "nested":
                        both = [(a, b, y) for a, b, y in zip(nested, intervals, truths)
                                if a is not None and b is not None]
                        paired.append(dict(
                            **identity, support=outcomes([a is not None for a in nested],
                                                        [b is not None for b in intervals]),
                            coverage=outcomes([covers(a, y) for a, _, y in both],
                                              [covers(b, y) for _, b, y in both]),
                            positive=outcomes([a[0] > 0 for a, _, _ in both],
                                              [b[0] > 0 for _, b, _ in both]),
                        ))
    return summary, paired


def write_table(output, name, rows, *, check):
    with (REFERENCE / name).open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields, expected = reader.fieldnames, list(reader)
    encoded = [
        {k: json.dumps(v) if isinstance(v, dict) else "" if v is None else str(v)
         for k, v in row.items()}
        for row in rows
    ]
    with (output / name).open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(encoded)
    if check:
        if len(encoded) != len(expected):
            raise ValueError(f"{name} has {len(encoded)} rows; expected {len(expected)}")
        for actual, saved in zip(encoded, expected):
            if saved.get("mean_width") and actual.get("mean_width"):
                if isclose(float(saved["mean_width"]), float(actual["mean_width"]),
                           rel_tol=1e-9, abs_tol=1e-12):
                    actual["mean_width"] = saved["mean_width"]
            if actual != saved:
                raise ValueError(f"{name} differs: {actual} != {saved}")


def run(output, workers, smoke):
    output.mkdir()
    plan = schedule(1 if smoke else 400)
    with (output / "intervals.jsonl").open("x", buffering=1) as stream:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = []
            for index, views in enumerate(pool.map(regenerate, plan), 1):
                for view in views:
                    stream.write(json.dumps(view, allow_nan=False) + "\n")
                    rows.append(view)
                if index % 100 == 0:
                    print(f"{index} native fields complete", flush=True)
    summary, paired = summarize(rows)
    if not smoke:
        if len(rows) != 2420:
            raise ValueError("Expected 2,400 native fields and 20 paired merge controls")
        write_table(output, "ALL_COMPARISONS.csv", summary, check=True)
        write_table(output, "PAIRED_COMPARISONS.csv", paired, check=True)
    print(f"{len(rows)} views regenerated" + (" (smoke run)" if smoke else "; both tables verified"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new output directory")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--smoke", action="store_true", help="first native field per condition and paired merges")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    run(args.output, args.workers, args.smoke)
