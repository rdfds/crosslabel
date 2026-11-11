"""Check adjacent tile sizes on the existing evaluation fields."""

import argparse
import json
from pathlib import Path

from field_only import DATA, block_lower, previous, rate


def run(field_results):
    saved = json.loads(field_results.read_text())
    original = previous.original
    rows = {row["seed"]: row for row in saved["rows"]
            if row["setting"] in original.SETTINGS}
    if len(rows) != 540:
        raise ValueError("Unexpected evaluation roster")
    result = {}
    for family, (setting, config) in enumerate(original.SETTINGS.items()):
        groups = {}
        for _, seed, group, replicate in previous.schedule(family, 2):
            row = rows.pop(seed)
            q = row["q"]
            fractions = original.ALTERNATIVE if q is None else original.prior.restricted(q)
            xy, codes, target = original.prior.generate(
                config, fractions, original.rng(family, 2, group, replicate)
            )
            if (setting, target) != (row["setting"], row["target"]):
                raise ValueError(f"Field mismatch for {seed}")
            lowers = {str(side): block_lower(xy, codes, side) for side in (4, 5, 6)}
            if not previous.equivalent(lowers["5"], row["dtc_lower"]):
                raise ValueError(f"Original tile score changed for {seed}")
            name = "alternative" if q is None else "adequate"
            groups.setdefault(name, []).append({
                "seed": seed, "q": q, "tile_lowers": lowers,
            })
        result[setting] = {
            name: {
                "decisions": {
                    side: rate([None if lower is None else lower > 0 for lower in
                                (row["tile_lowers"][side] for row in fields)])
                    for side in ("4", "5", "6")
                },
                "fields": fields,
            } for name, fields in groups.items()
        }
        print(f"Verified {setting}", flush=True)
    if rows:
        raise ValueError("Unvisited fields")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--field-results", type=Path, default=DATA / "field_only/FIELD_RESULTS.json")
    args = parser.parse_args()
    result = run(args.field_results)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
