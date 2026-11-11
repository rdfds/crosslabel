"""Regenerate and check the paired Crosslabel and Picasso SPINNA comparisons."""

import argparse
from collections import defaultdict
import csv
import json
from math import isclose
from pathlib import Path

import numpy as np

from crosslabel import field_check
from interval_fields import sparse_field
from replay_intervals import evaluate
import spinna_design
import spinna_library_design


DATA = Path(__file__).resolve().parent / "data"
FOLDERS = {"geometry": DATA / "spinna", "library": DATA / "spinna_library"}


def plans(study, paired_only=True):
    if study == "geometry":
        return {p["seed"]: p for p in spinna_design.schedule() if p["paired"] or not paired_only}
    return {p["seed"]: p for p in spinna_library_design.schedule()}


def inputs(study, plan):
    if study == "geometry":
        xy, codes, truth = spinna_design.generate(plan)
        return xy, codes, truth["m2"], spinna_design.FAMILIES[plan["family"]]["K"]
    xy, codes, truth = sparse_field(np.random.default_rng(plan["seed"]), plan["fractions"])
    return xy, codes, truth, 4


def replay_intervals(study, output, seed):
    schedule = plans(study, paired_only=False)
    if seed is not None and seed not in schedule:
        raise ValueError(f"{seed} is outside the {study} interval study")
    saved = defaultdict(dict)
    with (FOLDERS[study] / "ALL_INTERVALS.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if "/" in row["method"]:
                saved[int(row["key"])][row["mask"], row["method"]] = row
    if set(saved) != set(schedule):
        raise ValueError("Saved and regenerated interval rosters differ")
    output.mkdir()
    with (output / "intervals.jsonl").open("x", buffering=1) as stream:
        for plan in ([schedule[seed]] if seed is not None else schedule.values()):
            field_seed = plan["seed"]
            xy, codes, truth, code_count = inputs(study, plan)
            intervals = evaluate(xy, codes, 30., code_count, field_seed)
            for mask, methods in intervals.items():
                for mode, variants in methods.items():
                    for method, interval in variants.items():
                        key = mask, f"{mode}/{method}"
                        expected = saved[field_seed].pop(key)
                        if truth != float(expected["target"]):
                            raise ValueError(f"Target differs for {field_seed}")
                        if interval is None:
                            matches = expected["decision"] == "unsupported"
                        else:
                            positive = "reject_joint_zero_model" if interval[0] > 0 else "not_rejected"
                            covered = interval[0] <= truth <= interval[1]
                            matches = (
                                positive == expected["decision"] and str(covered) == expected["covered"]
                                and all(isclose(x, float(expected[k]), rel_tol=1e-9, abs_tol=1e-12)
                                        for x, k in zip(interval, ("lower", "upper")))
                            )
                        if not matches:
                            raise ValueError(f"Interval differs for {field_seed}, {key}")
            if saved[field_seed]:
                raise ValueError(f"Missing intervals for {field_seed}")
            stream.write(json.dumps(dict(seed=field_seed, target=truth, intervals=intervals)) + "\n")
    print(f"{1 if seed is not None else len(schedule)} interval fields verified")


def replay(study, output, seed):
    if study == "geometry":
        from spinna_comparison import run
    else:
        from spinna_library_comparison import main as run
    schedule = plans(study)
    if seed is not None and seed not in schedule:
        raise ValueError(f"{seed} is outside the {study} paired comparison")
    with (FOLDERS[study] / "PAIRED.csv").open(newline="") as stream:
        saved = {int(row["seed"]): row for row in csv.DictReader(stream)}
    if set(saved) != set(schedule):
        raise ValueError("Saved and regenerated comparison rosters differ")
    output.mkdir()
    for plan in ([schedule[seed]] if seed is not None else schedule.values()):
        field_seed = plan["seed"]
        xy, codes, truth, code_count = inputs(study, plan)
        case = output / str(field_seed)
        case.mkdir()
        source = case / "input.npz"
        np.savez_compressed(source, xy=xy, labels=codes)
        field = field_check(xy, codes, code_count=code_count, radius_nm=30., seed=field_seed)
        with (case / "crosslabel.json").open("x") as stream:
            json.dump(dict(target=truth, **field), stream, indent=2, allow_nan=False)
        decision = ("unsupported" if field["positive_lower_limit"] is None else
                    "reject_joint_zero_model" if field["positive_lower_limit"] else "not_rejected")
        expected = saved[field_seed]
        if decision != expected["nested"]:
            raise ValueError(f"Crosslabel decision differs for field {field_seed}")
        fits = case / "spinna"
        if study == "geometry":
            run(source, plan["family"], field_seed, fits)
            keys = ("approximate_omnibus_p", "approximate_improvement_p",
                    "reject_omnibus", "reject_improvement")
            expected_values = (float(expected["p_omnibus"]), float(expected["p_improvement"]),
                               expected["omnibus"] == "True", expected["improvement"] == "True")
        else:
            run(source, fits, field_seed)
            keys = ("model_bootstrap_p", "decision")
            expected_values = (float(expected["p"]), expected["spinna"])
        actual = json.loads((fits / "RESULT.json").read_text())
        if tuple(actual[k] for k in keys) != expected_values:
            raise ValueError(f"SPINNA ranks or decisions differ for field {field_seed}")
        print(f"{field_seed}: ranks and decisions match", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", choices=["geometry", "library"])
    parser.add_argument("output", type=Path, help="new output directory")
    parser.add_argument("--seed", type=int, help="one field from the study roster")
    parser.add_argument("--intervals-only", action="store_true",
                        help="verify all interval variants without SPINNA")
    args = parser.parse_args()
    command = replay_intervals if args.intervals_only else replay
    command(args.study, args.output, args.seed)
