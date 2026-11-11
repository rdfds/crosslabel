"""Regenerate the robustness fields and verify both methods' scores and decisions."""

import argparse
import hashlib
import json
from math import isclose
from pathlib import Path
import time

import numpy as np

from crosslabel import field_check

import robustness_design as original
from triple_correlation import score


ORIGINAL = Path(__file__).resolve().parent / "data/robustness_reference"
COMPARISON = ORIGINAL.parent / "robustness"

SCORES = ("dtc_connected", "dtc_connected_per_anchor", "dtc_raw",
          "crosslabel_point", "crosslabel_lower", "marked_pair")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def schedule(family, phase):
    compositions = original.GRID if phase == 1 else (*original.GRID, *original.OFF_GRID, None)
    for j, q in enumerate(compositions):
        size = 59 if phase == 1 else (80 if q is None else 20)
        for i in range(size):
            yield q, (phase + 10) * 100000 + 10000 * family + 100 * j + i, j, i


def equivalent(actual, saved):
    if isinstance(actual, dict) and isinstance(saved, dict):
        return actual.keys() == saved.keys() and all(equivalent(v, saved[k]) for k, v in actual.items())
    if isinstance(actual, (list, tuple)) and isinstance(saved, (list, tuple)):
        return len(actual) == len(saved) and all(equivalent(a, b) for a, b in zip(actual, saved))
    if isinstance(actual, float) and isinstance(saved, float):
        return isclose(actual, saved, rel_tol=1e-9, abs_tol=1e-12)
    return actual == saved


def train(family, settings, saved_training):
    endpoints = []
    for j, name in enumerate(("monomer", "dimer")):
        bank = []
        for i in range(64):
            xy, codes, _ = original.prior.generate(
                settings, original.prior.restricted(float(j)), original.rng(family, 0, j, i),
            )
            bank.append(original.summary(xy, codes, settings["K"]))
        endpoint = np.stack(bank)
        if not np.array_equal(endpoint, saved_training[name]):
            raise ValueError(f"Regenerated pair-training histogram differs for {name}, family {family}")
        endpoints.append(endpoint)
    return original.Library.train(*endpoints)


def collect_rows(setting, family, model, phase):
    folder = ORIGINAL / setting
    saved = original.records(folder / ("REFERENCE.jsonl" if phase == 1 else "EVALUATION.jsonl"))
    roster = list(schedule(family, phase))
    if len(saved) != len(roster):
        raise ValueError(f"Wrong saved roster length in {setting}")
    settings = original.SETTINGS[setting]
    for (q, seed, j, i), old in zip(roster, saved):
        if (q, seed) != (old["q"], old["seed"]):
            raise ValueError(f"Saved roster mismatch at {setting}, {seed}")
        fractions = original.ALTERNATIVE if q is None else original.prior.restricted(q)
        xy, codes, target = original.prior.generate(
            settings, fractions, original.rng(family, phase, j, i)
        )
        pair_q, pair_score = model.fit(original.summary(xy, codes, settings["K"]))
        if not equivalent((len(xy), float(target), float(pair_q[0]), float(pair_score[0])), (
                old["records"], float(old["target"]), old["pair_q"], old["pair_score"])):
            raise ValueError(f"Regenerated field parity failed at {setting}, {seed}")
        measured = field_check(xy, codes, code_count=settings["K"], radius_nm=30., seed=seed)
        interval = measured["interval_95_model"]
        measured_values = (
            measured["status"], measured["reason"], measured["estimate"], interval,
            None if interval is None else interval[0] <= target <= interval[1],
        )
        if not equivalent(measured_values, tuple(old[k] for k in
                                                 ("status", "reason", "estimate", "interval", "covered"))):
            raise ValueError(f"Regenerated Crosslabel interval differs at {setting}, {seed}")
        row = dict(setting=setting, q=q, seed=seed, records=len(xy), target=target,
                   covered=measured_values[-1], crosslabel_point=measured["estimate"],
                   crosslabel_lower=None if interval is None else interval[0], interval=interval,
                   marked_pair=float(pair_score[0]))
        try:
            row.update(score(xy, codes))
        except ValueError as exc:
            row.update(dtc_connected=None, dtc_connected_per_anchor=None,
                       dtc_raw=None, dtc_reason=str(exc))
        if phase == 2:
            row.update(old_lower_rank=old["nested_p"], old_pair_rank=old["pair_p"])
        yield row


def ranks(row, reference):
    values = {}
    for name in SCORES:
        bank = [[r[name] for r in reference if r["q"] == q] for q in original.GRID]
        if any(len(block) != 59 for block in bank):
            raise ValueError("Incomplete reference grid")
        values[name] = (None if row[name] is None or any(None in block for block in bank)
                        else float(max(original.rank_grid(row[name], bank))))
    if values["crosslabel_lower"] != row["old_lower_rank"] or values["marked_pair"] != row["old_pair_rank"]:
        raise ValueError(f"Prior ranks failed verification at {row['seed']}")
    return values


def collect(output):
    hashes = json.loads((ORIGINAL / "SHA256.json").read_text())
    for setting in original.SETTINGS:
        for filename in ("TRAINING.npz", "REFERENCE.jsonl", "EVALUATION.jsonl"):
            path = f"{setting}/{filename}"
            if sha(ORIGINAL / path) != hashes[path]:
                raise ValueError(f"Input checksum differs: {path}")
    output.mkdir(exist_ok=False)
    start = time.monotonic()
    expected = {
        name: original.records(COMPARISON / (name.upper() + ".jsonl"))
        for name in ("reference", "evaluation")
    }
    state = dict(status="running", reference=0, evaluation=0, completed=[])
    original.write(output / "STATE.json", state)
    try:
        with (output / "REFERENCE.jsonl").open("x", buffering=1) as ref, \
                (output / "EVALUATION.jsonl").open("x", buffering=1) as test:
            for family, setting in enumerate(original.SETTINGS):
                with np.load(ORIGINAL / setting / "TRAINING.npz", allow_pickle=False) as training:
                    model = train(family, original.SETTINGS[setting], training)
                reference = []
                for row in collect_rows(setting, family, model, 1):
                    if not equivalent(row, expected["reference"][state["reference"]]):
                        raise ValueError(f"Reference score differs at {setting}, {row['seed']}")
                    ref.write(json.dumps(row, allow_nan=False) + "\n")
                    reference.append(row)
                    state["reference"] += 1
                for row in collect_rows(setting, family, model, 2):
                    row["ranks"] = ranks(row, reference)
                    if not equivalent(row, expected["evaluation"][state["evaluation"]]):
                        raise ValueError(f"Evaluation score or rank differs at {setting}, {row['seed']}")
                    test.write(json.dumps(row, allow_nan=False) + "\n")
                    state["evaluation"] += 1
                state["completed"].append(setting)
                original.write(output / "STATE.json", state)
                print(f"{setting}: 177 reference and 180 evaluation fields scored", flush=True)
        if any(state[name] != len(rows) for name, rows in expected.items()):
            raise ValueError("Regenerated roster size differs")
        state["status"] = "complete"
    except BaseException as exc:
        state.update(status="failed_or_partial", error=repr(exc))
        raise
    finally:
        state["seconds"] = time.monotonic() - start
        original.write(output / "STATE.json", state)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    collect(parser.parse_args().output)
