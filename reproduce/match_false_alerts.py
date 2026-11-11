"""Describe operating points from the saved scores without refitting."""

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path


SCORES = (
    "dtc_connected",
    "dtc_connected_per_anchor",
    "crosslabel_point",
    "crosslabel_lower",
)
ON_GRID = (0.25, 0.50, 0.75)
OFF_GRID = (0.375, 0.625)
GROUPS = {q: 20 for q in (*ON_GRID, *OFF_GRID)} | {None: 80}
SETTINGS = (
    "balanced_lower_density",
    "low_retention_skewed_codes",
    "varying_intensity_skewed_codes",
)


def operating_points(rows, score, ordering):
    groups = defaultdict(Counter)
    for row in rows:
        value = -row[score] if ordering == "score" else row["ranks"][score]
        if not math.isfinite(value):
            raise ValueError(f"Nonfinite {score} value for seed {row['seed']}")
        groups[value][row["q"]] += 1

    counts = Counter()
    yield 0, 0, 0, 0
    for value in sorted(groups):
        counts.update(groups[value])
        on_grid = sum(counts[q] for q in ON_GRID)
        off_grid = sum(counts[q] for q in OFF_GRID)
        yield on_grid + off_grid, counts[None], on_grid, off_grid


def main(folder):
    source = (folder / "EVALUATION.jsonl").read_bytes()
    summary = json.loads((folder / "SUMMARY.json").read_text())
    if hashlib.sha256(source).hexdigest() != summary["sha256"]["EVALUATION.jsonl"]:
        raise ValueError("Saved scores do not match the comparison summary")
    rows = [json.loads(line) for line in source.splitlines()]
    if len(rows) != 540 or len({row["seed"] for row in rows}) != len(rows):
        raise ValueError("Expected 540 distinct evaluation fields")

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow((
        "setting", "method", "ordering", "false_alert_budget",
        "false_alerts", "on_grid_false_alerts", "off_grid_false_alerts",
        "detections",
    ))
    for setting in SETTINGS:
        subset = [row for row in rows if row["setting"] == setting]
        if Counter(row["q"] for row in subset) != GROUPS:
            raise ValueError(f"Unexpected composition roster for {setting}")
        for score in SCORES:
            for ordering in ("score", "reference_rank"):
                points = list(operating_points(subset, score, ordering))
                for budget in range(11):
                    false, detections, on_grid, off_grid = max(
                        point for point in points if point[0] <= budget
                    )
                    writer.writerow((
                        setting, score, ordering, budget,
                        false, on_grid, off_grid, detections,
                    ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, help="Saved comparison output directory")
    main(parser.parse_args().folder)
