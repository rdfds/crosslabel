"""Verify saved counts and regenerate the complete fixed sensitivity roster."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import binomtest
import retention_study as study


def equivalent(actual, expected):
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            equivalent(actual[k], v) for k, v in expected.items()
        )
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            equivalent(a, b) for a, b in zip(actual, expected)
        )
    if isinstance(expected, bool) or expected is None or isinstance(expected, str):
        return actual == expected
    if isinstance(expected, int):
        return actual == expected
    return math.isclose(actual, expected, rel_tol=1e-10, abs_tol=2e-12)


def main(replay):
    folder = Path(__file__).resolve().parent / "data/retention"
    rows = [json.loads(s) for s in (folder / "FIELDS.jsonl").read_text().splitlines()]
    summary = json.loads((folder / "SUMMARY.json").read_text())
    plans = list(study.schedule())
    if len(rows) != len(plans):
        raise ValueError("Incomplete sensitivity matrix")
    replays = []
    for row, plan in zip(rows, plans):
        if any(row[k] != v for k, v in plan.items()):
            raise ValueError("Changed seed or study setting")
        xy, codes, target = study.generate(plan)
        if study.catalogue_hash(xy, codes) != row["catalogue_sha256"] or target != row["target"]:
            raise ValueError("Regenerated catalogue or target differs")
        interval = row["result"]["interval_95_model"]
        if interval is None or row["flag"] != (interval[0] > 0):
            raise ValueError("Saved interval decision differs")
        if replay and plan["replicate"] == 1:
            result = study.field_check(xy, codes, code_count=4, radius_nm=30., seed=plan["seed"])
            if not equivalent(result, row["result"]) or result["positive_lower_limit"] != row["flag"]:
                raise ValueError("Replayed interval or decision differs")
            replays.append(plan["seed"])
    for cell in summary["cells"]:
        group = [r for r in rows if (r["retention"], r["expected_records"], r["alternative"])
                 == (cell["retention"], cell["expected_records"], cell["alternative"])]
        k = sum(r["flag"] for r in group)
        ci = list(binomtest(k, len(group)).proportion_ci(method="exact"))
        if len(group) != 80 or k != cell["alerts"] or not equivalent(ci, cell["ci95"]):
            raise ValueError("Saved alert count or exact limits differ")
        if sum(r["covered"] for r in group) != cell["covered"]:
            raise ValueError("Saved coverage differs")
    print(json.dumps(dict(catalogues_regenerated=len(rows), cell_tables_checked=len(summary["cells"]),
                          field_intervals_replayed=len(replays), status="pass"), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", action="store_true")
    main(parser.parse_args().replay)
