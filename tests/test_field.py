import numpy as np
import pytest

from crosslabel.field import field_check
from crosslabel.field_screen import density_score, reference_tail


def test_density_screen_flags_concentration_at_fixed_reference_tail():
    x = np.arange(8, dtype=float)
    regular = np.column_stack([np.repeat(x, 8), np.tile(x, 8)])
    concentrated = np.vstack([regular, np.repeat([[0., 0.]], 64, axis=0)])
    assert density_score(regular)["score"] == 0.
    assert reference_tail(density_score(concentrated)["score"], np.zeros(99)) == .01
    assert reference_tail(density_score(regular)["score"], np.zeros(99)) == 1.


def test_field_entry_point_keeps_unsupported_inputs_distinct_from_negative():
    xy = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    with pytest.raises(ValueError, match="two-record minimum"):
        field_check(xy, np.arange(4), code_count=4, radius_nm=30, seed=1)
