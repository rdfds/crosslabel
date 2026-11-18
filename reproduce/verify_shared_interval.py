"""Verify all saved shared-interval points and replay a fixed subset."""
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from crosslabel import _field_kernel as kernel
from collect_robustness import equivalent
from dtc_interval import spatial, value_gradient

folder = Path(__file__).resolve().parent / "data/shared_interval"
rows = [json.loads(s) for s in (folder / "FIELDS.jsonl").read_text().splitlines()]
if len(rows) != 540:
    raise ValueError("Complete 540-field roster required")
replayed = []
for row in rows:
    with np.load(folder / "arrays" / f'{row["seed"]}.npz') as arrays:
        objective = lambda x, W: value_gradient(x, W, float(arrays["area"]), float(arrays["disc"]))
        point, _ = objective(arrays["fine"].sum(0), arrays["W"])
        if not equivalent(float(point), row["result"]["value"]):
            raise ValueError("Saved point differs")
        if row["seed"] % 10000 not in (*range(5), *range(500, 505)):
            continue
        scale = spatial(arrays["fine"], arrays["W"], tuple(arrays["shape"]),
                        int(arrays["coarse_blocks"]), objective)
        with patch.object(kernel, "value_gradient", objective):
            result, calibration = kernel.nested_estimate(
                arrays["coarse"], arrays["W"], scale, row["seed"])
        if not equivalent(result, row["result"]):
            raise ValueError("Shared interval replay differs")
        for name, values in calibration.items():
            np.testing.assert_allclose(values, arrays[name], rtol=1e-10, atol=2e-12)
        replayed.append(row["seed"])
if len(replayed) != 30:
    raise ValueError("Wrong replay roster")
print("All 540 point scores and 30 full nested interval replays verified")
