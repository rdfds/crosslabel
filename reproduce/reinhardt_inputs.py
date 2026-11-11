"""Verify named Figure 2 arrays or convert extracted HDF5 members offline."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PIXEL_NM = 130.0


def fingerprint(array):
    array = np.ascontiguousarray(array)
    return dict(dtype=array.dtype.str, shape=list(array.shape), sha256=hashlib.sha256(array.tobytes()).hexdigest())


def convert(path, kind, rows):
    import h5py

    with h5py.File(path, "r") as f:
        if not isinstance(f.get("locs", getlink=True), h5py.HardLink):
            raise ValueError(f"Expected a local locs table: {path}")
        d = f["locs"]
        if not isinstance(d, h5py.Dataset) or d.is_virtual or d.external or d.shape != (rows,):
            raise ValueError(f"Unexpected table geometry: {path}")
        a = d[:]
    lpx, lpy = a["lpx"], a["lpy"]
    if kind == "locs":
        lpx, lpy = lpx.astype(float), lpy.astype(float)
    out = dict(x=(a["x"]*PIXEL_NM).astype("f4"), y=(a["y"]*PIXEL_NM).astype("f4"),
               lp=(PIXEL_NM * np.sqrt((lpx**2 + lpy**2)/2)).astype("f4"))
    if kind == "locs":
        out.update(frame=a["frame"].astype("i4"), group=a["group"].astype("i4"))
    else:
        out["n"] = a["n"]
    return out


def expected_files():
    return json.loads((HERE / "NUP96_ARRAYS.json").read_text())


def verify_npz(cache):
    files = {}
    for name, expected in expected_files().items():
        path = cache / name
        with np.load(path, allow_pickle=False) as data:
            if {k: fingerprint(data[k]) for k in data.files} != expected:
                raise ValueError(f"Reduced array identity differs: {name}")
        files[name] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size)
    return dict(cache=str(cache.resolve()), files=files)


def prepare(original, output):
    output.mkdir()
    for name, expected in expected_files().items():
        _, kind, round_file = name.split("-")
        round_id = round_file.split(".")[0]
        if kind == "locs":
            path = original / "RESI" / "RESI rounds" / f"RESI_{round_id}_apicked.hdf5"
        else:
            path = original / "RESI" / "RESI result" / f"RESI_{round_id}_apicked_resi.hdf5"
        arrays = convert(path, kind, expected["x"]["shape"][0])
        if {k: fingerprint(v) for k, v in arrays.items()} != expected:
            raise ValueError(f"Extracted member differs from the retained data: {name}")
        np.savez_compressed(output / name, **arrays)
    (output / "MANIFEST.json").write_text(json.dumps(verify_npz(output), indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figure2", type=Path, help="Extracted Figure2,EDF06,EDF07 directory")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    prepare(args.figure2, args.output)
