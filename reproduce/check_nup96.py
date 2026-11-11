"""Verify the frozen Nup96 comparison and summarize its paired decisions."""

import argparse
import hashlib
import json
from math import isclose
from pathlib import Path

from scipy.stats import binomtest
from crosslabel.field_screen import reference_tail

SCORES = ("dtc_raw", "dtc_connected", "crosslabel_point", "crosslabel_lower")


def jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def proportion(positive, total):
    ci = binomtest(positive, total).proportion_ci() if total else None
    return {"positive": positive, "supported": total,
            "exact_95_ci": None if ci is None else [ci.low, ci.high]}


def check_summary(saved, result):
    for method, summary in result["methods"].items():
        for cohort in ("two_code", "four_code"):
            previous = saved["methods"][method][cohort]
            old = previous["exact_95_ci"]
            new = summary[cohort]["exact_95_ci"]
            if (isinstance(old, list) and new is not None and len(old) == len(new)
                    and all(isclose(a, b, rel_tol=1e-9) for a, b in zip(old, new))):
                previous["exact_95_ci"] = new
    if saved != result:
        raise ValueError("Recomputed summary differs from the saved result")


def summarize(output, check=False):
    state = json.loads((output / "STATE.json").read_text())
    refs = jsonl(output / "REFERENCE.jsonl")
    rows = jsonl(output / "PAIRED.jsonl")
    if state["status"] != "complete" or (len(refs), len(rows)) != (99, 200):
        raise ValueError("The complete frozen matrix is required")
    if [r["seed"] for r in refs] != list(range(988000, 988099)):
        raise ValueError("Reference roster mismatch")
    expected = [(seed, adequate) for seed in range(987000, 987100)
                for adequate in (True, False)]
    if [(r["seed"], r["adequate"]) for r in rows] != expected:
        raise ValueError("Evaluation roster mismatch")

    summary = {}
    for name in SCORES:
        bank = [r[name] for r in refs]
        cohorts = {}
        for adequate, label in ((True, "two_code"), (False, "four_code")):
            selected = [r for r in rows if r["adequate"] == adequate]
            supported = []
            for r in selected:
                rank = r["ranks"][name]
                expected_rank = None if r[name] is None else reference_tail(r[name], bank)
                if rank != expected_rank:
                    raise ValueError("Saved rank differs from the frozen rule")
                if rank is not None:
                    supported.append(rank <= .05)
            paired = [(r["ranks"][name] <= .05, r["ranks"]["crosslabel_lower"] <= .05)
                      for r in selected if r["ranks"][name] is not None
                      and r["ranks"]["crosslabel_lower"] is not None]
            only_method = paired.count((True, False))
            only_crosslabel = paired.count((False, True))
            discordant = only_method + only_crosslabel
            cohorts[label] = dict(
                planned=len(selected), unsupported=len(selected) - len(supported),
                **proportion(sum(supported), len(supported)),
                compared_with_crosslabel_lower=len(paired),
                both_positive=paired.count((True, True)),
                both_negative=paired.count((False, False)),
                only_method=only_method, only_crosslabel_lower=only_crosslabel,
                exact_paired_p=(binomtest(only_method, discordant).pvalue if discordant else 1.0),
            )
        score_pairs = [(a[name], b[name]) for a, b in zip(rows[::2], rows[1::2])
                       if a[name] is not None and b[name] is not None]
        summary[name] = {
            "reference_count": len(bank), "reference_min": min(bank), "reference_max": max(bank),
            "alternative_higher": sum(b > a for a, b in score_pairs),
            "control_higher": sum(a > b for a, b in score_pairs),
            "ties": sum(a == b for a, b in score_pairs),
            **cohorts,
        }
    paths = ("REFERENCE.jsonl", "PAIRED.jsonl", "NATIVE.json", "STATE.json")
    result = {
        "native": json.loads((output / "NATIVE.json").read_text()),
        "methods": summary,
        "collection": state,
        "sha256": {p: hashlib.sha256((output / p).read_bytes()).hexdigest() for p in paths},
    }
    summary_path = output / "SUMMARY.json"
    if check:
        check_summary(json.loads(summary_path.read_text()), result)
    else:
        with summary_path.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    for name, s in summary.items():
        print(name, "control alerts", s["two_code"]["positive"],
              "alternative alerts", s["four_code"]["positive"],
              "alt score higher", s["alternative_higher"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    summarize(args.output, check=args.check)
