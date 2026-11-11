"""Observation settings and random streams for the robustness comparison."""

import json

import numpy as np

import simulation as prior
from marked_pairs import Library, rank_grid, summary


BASE = dict(shape="iid", spacing=7., K=4, mask=False)


SETTINGS = {
    "balanced_lower_density": dict(BASE, density=38., p=.60, codes=(.25, .25, .25, .25)),
    "low_retention_skewed_codes": dict(BASE, density=55., p=.27, codes=(.10, .20, .30, .40)),
    "varying_intensity_skewed_codes": dict(BASE, density=65., p=.47, codes=(.15, .26, .48, .11),
                                           mask=True, mask_pattern="two_dimensional"),
}


GRID = (.25, .50, .75)


OFF_GRID = (.375, .625)


ALTERNATIVE = (.65, .20, .15)


ALPHA = .05


def write(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def rng(family, phase, composition, repeat):
    return np.random.default_rng(np.random.SeedSequence([461902, family, phase, composition, repeat]))


def records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]
