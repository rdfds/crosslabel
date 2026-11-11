"""Run the documented field commands against the packaged example."""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from crosslabel.field_cli import main


EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "field"


def check_report(output):
    expected = json.loads((EXAMPLES / "matched_alternative.json").read_text())
    actual = json.loads(output.read_text())
    for name in ("status", "positive_lower_limit", "support", "seed", "radius_nm"):
        assert actual[name] == expected[name]
    for name in ("estimate", "interval_95_model", "window_nm"):
        np.testing.assert_allclose(actual[name], expected[name], rtol=1e-10, atol=1e-10)


def test_csv_example(tmp_path):
    output = tmp_path / "report.json"
    main(["--csv", str(EXAMPLES / "matched_alternative.csv"),
          "--radius-nm", "30", "--seed", "984500", "--output", str(output)])
    check_report(output)


def test_picasso_example(tmp_path):
    h5py = pytest.importorskip("h5py")
    with (EXAMPLES / "matched_alternative.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    paths = []
    for code in sorted({row["code"] for row in rows}):
        path = tmp_path / f"round-{code}.hdf5"
        selected = [(float(row["x_nm"]), float(row["y_nm"]))
                    for row in rows if row["code"] == code]
        with h5py.File(path, "w") as handle:
            handle["locs"] = np.array(selected, dtype=[("x", "f8"), ("y", "f8")])
        paths.append(str(path))
    output = tmp_path / "report.json"
    main(["--picasso", *paths, "--pixel-nm", "1", "--radius-nm", "30",
          "--seed", "984500", "--output", str(output)])
    check_report(output)
