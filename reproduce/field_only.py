"""Evaluate fixed field-only decisions and a controlled intensity perturbation."""

import argparse
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
from scipy.stats import binomtest, t

from crosslabel._field_kernel import fine_sufficient, value_gradient
from crosslabel.field import field_check

import collect_robustness as previous


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OLD = DATA / "robustness"
CELLS = {"OBZ_2D": "cell1_2023_Pos3", "OBZ_3D": "cell4_3D"}

def block_lower(xy, codes, side=5):
    lo = xy.min(axis=0)
    span = np.ptp(xy, axis=0)
    if np.any(span == 0):
        return None
    bins = np.minimum(((xy - lo) / span * side).astype(int), side - 1)
    tile = bins @ [side, 1]
    count = side**2
    try:
        scores = np.array([
            previous.score(xy[tile == j], codes[tile == j])["dtc_connected"]
            for j in range(count)
        ])
    except ValueError:
        return None
    return float(scores.mean() - t.ppf(.975, count - 1) * scores.std(ddof=1) / side)


def evaluate(robustness=OLD):
    original = previous.original
    saved = original.records(robustness / "EVALUATION.jsonl")
    if len(saved) != 540 or len({r["seed"] for r in saved}) != 540:
        raise ValueError("Unexpected evaluation roster")
    by_seed = {r["seed"]: r for r in saved}
    rows = []
    for family, (setting, config) in enumerate(original.SETTINGS.items()):
        for q, seed, group, replicate in previous.schedule(family, 2):
            row = by_seed[seed]
            if (row["setting"], row["q"]) != (setting, q):
                raise ValueError(f"Unexpected setting or composition for {seed}")
            fractions = original.ALTERNATIVE if q is None else original.prior.restricted(q)
            rng = lambda: original.rng(family, 2, group, replicate)
            xy, codes, target = original.prior.generate(config, fractions, rng())
            direct = previous.score(xy, codes)["dtc_connected"]
            if (len(xy), target) != (row["records"], row["target"]) or not np.isclose(
                direct, row["dtc_connected"], rtol=1e-12, atol=1e-12
            ):
                raise ValueError(f"Saved field mismatch for {seed}")
            rows.append(dict(
                setting=setting, q=q, seed=seed, target=target,
                point=row["crosslabel_point"], interval=row["interval"],
                dtc_lower=block_lower(xy, codes),
            ))
            if setting != "varying_intensity_skewed_codes":
                continue
            with patch.object(original.prior.old, "mask_pmf", return_value=np.ones((20, 20))):
                uniform, uniform_codes, uniform_target = original.prior.generate(
                    dict(config, mask_pattern="uniform_diagnostic"), fractions, rng()
                )
            if not np.array_equal(codes, uniform_codes) or uniform_target != target:
                raise ValueError(f"Control altered marks or target for {seed}")
            result = field_check(uniform, codes, code_count=4, radius_nm=30., seed=seed)
            rows.append(dict(
                setting="uniform_matched_to_varying", q=q, seed=seed, target=target,
                point=result["estimate"], interval=result["interval_95_model"],
                dtc_lower=block_lower(uniform, codes),
            ))
        print(f"Completed {setting}", flush=True)
    return rows


def rate(flags):
    supported = [flag for flag in flags if flag is not None]
    k, n = sum(supported), len(supported)
    ci = binomtest(k, n).proportion_ci() if n else (None, None)
    return dict(calls=k, supported=n, total=len(flags), ci95=list(ci))


def summarize(rows):
    result = {}
    for setting in dict.fromkeys(r["setting"] for r in rows):
        groups = {}
        for group, adequate in (("adequate", True), ("alternative", False)):
            selected = [r for r in rows if r["setting"] == setting and (r["q"] is not None) == adequate]
            intervals = [r for r in selected if r["interval"] is not None]
            groups[group] = dict(
                crosslabel=rate([None if r["interval"] is None else r["interval"][0] > 0 for r in selected]),
                dtc_tiles=rate([None if r["dtc_lower"] is None else r["dtc_lower"] > 0 for r in selected]),
                coverage=rate([
                    None if r["interval"] is None else r["interval"][0] <= r["target"] <= r["interval"][1]
                    for r in selected
                ]),
                mean_error=float(np.mean([r["point"] - r["target"] for r in intervals])),
                median_lower_offset=float(np.median([r["point"] - r["interval"][0] for r in intervals])),
                median_upper_offset=float(np.median([r["interval"][1] - r["point"] for r in intervals])),
                misses_below=sum(r["target"] < r["interval"][0] for r in intervals),
                misses_above=sum(r["target"] > r["interval"][1] for r in intervals),
            )
        result[setting] = groups
    return result


def verified(path, digest):
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"Input checksum mismatch: {path.name}")
    return path


def source_xy(path):
    with h5py.File(path, "r") as data:
        locs = data["locs"][:]
        return (np.column_stack((locs["x"], locs["y"])) * 130.).astype(float)


def g5m_xy(root, cell, number, manifest):
    name = f"g5m-{cell}-g5mmol-R{number}.npz"
    entry = manifest[name]
    if (root / name).is_file():
        with np.load(verified(root / name, entry["sha256"])) as data:
            return np.column_stack((data["x"], data["y"])).astype(float)
    path = root / "data/Fig4" / cell / f"R{number}_picked_dbscan_molmap.hdf5"
    xy = source_xy(path)
    if hashlib.sha256(xy.astype("<f8", copy=False).tobytes()).hexdigest() != entry["xy_sha256"]:
        raise ValueError(f"Scored coordinates differ: {path.name}")
    return xy


def catalogue_score(rounds):
    xy = np.concatenate(rounds)
    counts = [len(x) for x in rounds]
    codes = np.repeat(np.arange(4), counts)
    stats, _, _ = fine_sufficient(xy, codes, 30., 4)
    all_codes = stats["all"]
    point, _ = value_gradient(all_codes["x"].sum(axis=0), all_codes["W"])
    return dict(records=sum(counts), round_counts=counts, point=float(point)), xy, codes


def g5m(native, ga, sensitivity=False):
    reduction_manifest = json.loads((HERE / "CD20_INPUTS.json").read_text())
    ga_manifest = json.loads((DATA / "source_catalogues/GA_SHA256.json").read_text())
    result = {}
    for i, (name, cell) in enumerate(CELLS.items()):
        result[cell] = {}
        catalogs = {}
        for method in ("G5M", "GA"):
            rounds = []
            for r in range(1, 5):
                if method == "G5M":
                    xy = g5m_xy(native, cell, r, reduction_manifest)
                else:
                    file = f"{name}-R{r}.hdf5"
                    path = (ga / file if ga is not None else
                            native / "data/Fig4" / cell / f"R{r}_picked_ga_centers.hdf5")
                    xy = source_xy(verified(path, ga_manifest[file]))
                rounds.append(xy)
            catalogs[method] = rounds
            summary, xy, codes = catalogue_score(rounds)
            if sensitivity:
                summary["connected_dtc_per_anchor"] = previous.score(xy, codes)["dtc_connected_per_anchor"]
            result[cell][method] = summary
        if sensitivity:
            seed = 125831 + i
            rng = np.random.default_rng(seed)
            sizes = result[cell]["GA"]["round_counts"]
            draws = []
            for _ in range(32):
                rounds = [
                    xy[np.sort(rng.choice(len(xy), size=n, replace=False))]
                    for xy, n in zip(catalogs["G5M"], sizes)
                ]
                draws.append(catalogue_score(rounds)[0]["point"])
            result[cell]["count_matched_g5m"] = dict(seed=seed, round_counts=sizes, points=draws)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--cd20-only", action="store_true")
    mode.add_argument("--cd20-sensitivity", action="store_true")
    parser.add_argument("--g5m-dir", type=Path)
    parser.add_argument("--ga-dir", type=Path)
    parser.add_argument("--robustness-dir", type=Path, default=OLD,
                        help="directory with saved or regenerated robustness scores")
    args = parser.parse_args()
    if args.cd20_only or args.cd20_sensitivity:
        if args.g5m_dir is None:
            parser.error("the CD20 comparison requires --g5m-dir")
        result = {"obinutuzumab_cd20": g5m(args.g5m_dir, args.ga_dir, args.cd20_sensitivity)}
    else:
        rows = evaluate(args.robustness_dir)
        result = {"summary": summarize(rows), "rows": rows}
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
