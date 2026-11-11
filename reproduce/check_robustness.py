"""Verify all frozen-grid ranks and report detection and false alerts."""

import argparse
import csv
import json
from math import isclose
from pathlib import Path

from scipy.stats import binomtest

from collect_robustness import SCORES, original, ranks, schedule, sha


def check_rates(saved, recomputed):
    if len(saved) != len(recomputed):
        raise ValueError("Rate row counts differ")
    for number, (old, new) in enumerate(zip(saved, recomputed), start=1):
        if old.keys() != new.keys():
            raise ValueError(f"Rate columns differ in row {number}")
        for key, value in old.items():
            if value == new[key]:
                continue
            if (key in ("ci_low", "ci_high") and value not in (None, "")
                    and new[key] not in (None, "")
                    and isclose(float(value), float(new[key]), rel_tol=1e-9)):
                continue
            raise ValueError(f"Rate row {number} differs in {key}: {value!r} vs {new[key]!r}")


def rates(rows, score):
    decisions = [r["ranks"][score] <= original.ALPHA for r in rows
                 if r["ranks"][score] is not None]
    count = sum(decisions)
    interval = binomtest(count, len(decisions)).proportion_ci() if decisions else None
    paired = [(r["ranks"][score] <= original.ALPHA,
               r["ranks"]["crosslabel_lower"] <= original.ALPHA)
              for r in rows if r["ranks"][score] is not None
              and r["ranks"]["crosslabel_lower"] is not None]
    only_score = paired.count((True, False))
    only_lower = paired.count((False, True))
    discordant = only_score + only_lower
    return dict(planned=len(rows), supported=len(decisions),
                unsupported=len(rows) - len(decisions), alerts=count,
                ci_low=None if interval is None else interval.low,
                ci_high=None if interval is None else interval.high,
                paired_n=len(paired), both_alert=paired.count((True, True)),
                neither_alert=paired.count((False, False)),
                only_score=only_score, only_lower=only_lower,
                paired_p=binomtest(only_score, discordant).pvalue if discordant else 1.0)


def summarize(folder):
    state = json.loads((folder / "STATE.json").read_text())
    reference = original.records(folder / "REFERENCE.jsonl")
    evaluation = original.records(folder / "EVALUATION.jsonl")
    if state["status"] != "complete" or (len(reference), len(evaluation)) != (531, 540):
        raise ValueError("The complete fixed matrix is required")
    all_rows = []
    for family, setting in enumerate(original.SETTINGS):
        ref = [r for r in reference if r["setting"] == setting]
        test = [r for r in evaluation if r["setting"] == setting]
        for phase, rows in ((1, ref), (2, test)):
            expected = [(q, seed) for q, seed, _, _ in schedule(family, phase)]
            if [(r["q"], r["seed"]) for r in rows] != expected:
                raise ValueError(f"Roster mismatch in {setting}, phase {phase}")
        for row in test:
            if ranks(row, ref) != row["ranks"]:
                raise ValueError(f"Saved ranks differ at {row['seed']}")
        groups = [
            (str(q), [r for r in test if r["q"] == q])
            for q in (*original.GRID, *original.OFF_GRID)
        ]
        groups.extend([
            ("on_grid", [r for r in test if r["q"] in original.GRID]),
            ("off_grid", [r for r in test if r["q"] in original.OFF_GRID]),
            ("alternative", [r for r in test if r["q"] is None]),
        ])
        for group, rows in groups:
            for score in SCORES:
                all_rows.append(dict(setting=setting, group=group, score=score,
                                     **rates(rows, score)))
    return {
        "collection": state,
        "sha256": {name: sha(folder / name) for name in
                   ("REFERENCE.jsonl", "EVALUATION.jsonl", "STATE.json")},
        "rates": all_rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = summarize(args.output)
    path = args.output / "SUMMARY.json"
    if args.check:
        saved = json.loads(path.read_text())
        check_rates(saved.pop("rates"), result["rates"])
        if saved != {key: value for key, value in result.items() if key != "rates"}:
            raise ValueError("Recomputed collection state or input hashes differ")
        with (args.output / "RATES.csv").open() as handle:
            check_rates(list(csv.DictReader(handle)), [
                {k: "" if v is None else str(v) for k, v in row.items()}
                for row in result["rates"]
            ])
    else:
        with path.open("x") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write("\n")
        with (args.output / "RATES.csv").open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=result["rates"][0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(result["rates"])
    print("Verified 531 reference fields, 540 evaluation fields and 144 rate rows.")
