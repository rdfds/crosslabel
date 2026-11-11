"""Regenerate the Nup96 results from verified localization catalogues."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np

from crosslabel.field import field_check
from crosslabel.field_screen import reference_tail
import reinhardt_inputs as inputs
from triple_correlation import score


DATA = Path(__file__).resolve().parent / "data"
SAVED = DATA / "nup96"
SCORES = ("dtc_raw", "dtc_connected", "crosslabel_point", "crosslabel_lower")


def jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def row_for(xy, labels, crosslabel):
    lower = None if crosslabel["interval_95_model"] is None else crosslabel["interval_95_model"][0]
    row = dict(crosslabel_point=crosslabel["estimate"], crosslabel_lower=lower,
               crosslabel_interval=crosslabel["interval_95_model"], records=len(xy),
               code_counts=np.bincount(labels, minlength=4).tolist())
    try:
        row.update(score(xy, labels))
    except ValueError as exc:
        row.update(dtc_raw=None, dtc_connected=None, dtc_reason=str(exc))
    return row


def load_templates(directory):
    def load(name):
        return np.load(directory / name, allow_pickle=False)

    centres, centre_codes, pieces = [], [], []
    with (DATA / "nup96_groups/ROWS.csv").open() as f:
        roster = [int(r["group"]) for r in csv.DictReader(f) if r["matching_support"] == "True"]
    if len(roster) != 1206:
        raise ValueError("Unexpected retained group roster")
    for r in range(4):
        with load(f"npc-resi-R{r+1}.npz") as z:
            centres.append(np.column_stack([z["x"], z["y"]]))
            centre_codes.append(np.full(len(z["x"]), r, dtype=int))
        with load(f"npc-locs-R{r+1}.npz") as z:
            keep = z["lp"] <= 8
            xy = np.column_stack([z["x"][keep], z["y"][keep]]).astype(float)
            groups = z["group"][keep]
            order = np.argsort(groups, kind="stable")
            cuts = np.searchsorted(groups[order], np.arange(1219))
            pieces.append((xy[order], cuts))
    templates = []
    rng = np.random.default_rng(986500)
    for group in roster:
        parts, codes = [], []
        for r, (xy, cuts) in enumerate(pieces):
            points = xy[cuts[group]:cuts[group+1]]
            selected = rng.choice(len(points), size=min(8, len(points)), replace=False)
            parts.append(points[selected])
            codes.append(np.full(len(selected), r, dtype=int))
        origin = np.concatenate([xy[cuts[group]:cuts[group+1]] for xy, cuts in pieces]).mean(0)
        templates.append((np.concatenate(parts) - origin, np.concatenate(codes)))
    return np.concatenate(centres), np.concatenate(centre_codes), templates


def template_pair(templates, seed):
    rng = np.random.default_rng(seed)
    count = rng.poisson(180)
    centres = rng.uniform(0, 10000, (count, 2))
    indices = rng.integers(len(templates), size=count)
    pos, alt, null = [], [], []
    witnesses = 0
    for center, index in zip(centres, indices):
        xy, marks = templates[index]
        angle = rng.uniform(0, 2 * np.pi)
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        x = xy @ rotation.T + center
        a = rng.permutation(4)[marks]
        b = rng.choice(4, size=2, replace=False)[a // 2]
        close = np.sum((x[:, None, :] - x[None, :, :])**2, axis=-1) <= 30**2
        neighbours = close @ (a[:, None] == np.arange(4))
        neighbours[np.arange(len(a)), a] = 0
        witnesses += int(((neighbours.sum(1)**2 - (neighbours**2).sum(1)) // 2).sum())
        if np.unique(b).size > 2:
            raise ValueError("Constructed null contains more than two codes")
        pos.append(x)
        alt.append(a)
        null.append(b)
    return np.concatenate(pos), np.concatenate(alt), np.concatenate(null), count, witnesses


def reinhardt(source, arrays):
    if source is None:
        inputs.verify_npz(arrays)
        return load_templates(arrays)
    with TemporaryDirectory(prefix="crosslabel-nup96-") as folder:
        converted = Path(folder) / "arrays"
        inputs.prepare(source, converted)
        return load_templates(converted)


def fit(xy, codes, seed):
    return field_check(xy, codes, code_count=4, radius_nm=30., seed=seed)


def zaza_row(path):
    metadata = json.loads((DATA / "source_catalogues/ZAZA_DESCRIPTIVE.json").read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["source_sha256"]:
        raise ValueError("Zaza catalogue checksum differs")
    with h5py.File(path, "r") as data:
        locs = data["locs"]
        xy = np.column_stack([locs["x"], locs["y"]]).astype(float) * metadata["pixel_nm"]
        codes = locs["resi_channel_id"].astype(int)
    result = field_check(xy, codes, code_count=4, radius_nm=30., seed=metadata["seed"])
    return row_for(xy, codes, result)


def check(expected, actual, label):
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError(f"Fields differ for {label}")
        for key in expected:
            check(expected[key], actual[key], f"{label}/{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f"Lengths differ for {label}")
        for i, (old, new) in enumerate(zip(expected, actual)):
            check(old, new, f"{label}/{i}")
    elif isinstance(expected, float) and not label.startswith("ranks/"):
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10, err_msg=label)
    elif actual != expected:
        raise ValueError(f"Value differs for {label}: {actual!r} != {expected!r}")


def replay(source, arrays, zaza):
    xy, codes, templates = reinhardt(source, arrays)
    native = {
        "Reinhardt": row_for(xy, codes, fit(xy, codes, 985000)),
        "Zaza": zaza_row(zaza),
    }
    check(json.loads((SAVED / "NATIVE.json").read_text()), native, "native")
    print("Both native catalogues reproduced", flush=True)

    refs = jsonl(SAVED / "REFERENCE.jsonl")
    paired = jsonl(SAVED / "PAIRED.jsonl")
    if [row["seed"] for row in refs] != list(range(988000, 988099)):
        raise ValueError("Reference roster differs")
    if [(row["seed"], row["adequate"]) for row in paired] != [
        (seed, adequate) for seed in range(987000, 987100) for adequate in (True, False)
    ]:
        raise ValueError("Paired roster differs")

    bank = {name: [] for name in SCORES}
    for old in refs:
        seed = old["seed"]
        xy, _, control, centers, _ = template_pair(templates, seed)
        actual = dict(seed=seed, centers=centers,
                      **row_for(xy, control, fit(xy, control, seed)))
        check(old, actual, f"reference/{seed}")
        for name in SCORES:
            bank[name].append(actual[name])
    print("99 reference fields reproduced", flush=True)

    for first, second in zip(paired[::2], paired[1::2]):
        seed = first["seed"]
        xy, alternative, control, centers, witnesses = template_pair(templates, seed)
        for old, labels, witness_count in ((first, control, 0), (second, alternative, witnesses)):
            actual = dict(seed=seed, adequate=old["adequate"], centers=centers, witnesses=witness_count,
                          **row_for(xy, labels, fit(xy, labels, seed)))
            ranks = {name: reference_tail(actual[name], bank[name]) for name in SCORES}
            check(old["ranks"], ranks, "ranks")
            check({key: value for key, value in old.items() if key != "ranks"},
                  actual, f"paired/{seed}/{old['adequate']}")
    print("200 paired fields and all rank decisions reproduced")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--reinhardt", type=Path, help="Extracted Figure2,EDF06,EDF07 directory")
    source.add_argument("--reinhardt-arrays", type=Path, help="Previously reduced, checksum-matched arrays")
    parser.add_argument("--zaza", required=True, type=Path, help="Original Nup96 Picasso HDF5")
    args = parser.parse_args()
    replay(args.reinhardt, args.reinhardt_arrays, args.zaza)
