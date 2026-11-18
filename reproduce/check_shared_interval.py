"""Check score identity, derivatives and shared spatial covariance."""
import importlib.util
from pathlib import Path
import sys

import numpy as np

from crosslabel import _field_kernel as kernel
from dtc_interval import sufficient, spatial, value_gradient

import robustness_design as benchmark
from triple_correlation import score


def main():
    for family, (setting, config) in enumerate(benchmark.SETTINGS.items()):
        xy, codes, _ = benchmark.prior.generate(
            config, benchmark.ALTERNATIVE, benchmark.rng(family, 2, 5, 0)
        )
        x, mapping, shape, blocks, W, area, disc = sufficient(xy, codes)
        objective = lambda sums, weights: value_gradient(sums, weights, area, disc)
        point, gradient = objective(x.sum(0), W)
        np.testing.assert_allclose(point, score(xy, codes)["dtc_connected_per_anchor"],
                                   rtol=1e-10, atol=2e-12)
        sums = x.sum(0)
        numerical = []
        for i in range(len(sums)):
            step = 1e-5*max(abs(sums[i]), .001)
            plus, minus = sums.copy(), sums.copy()
            plus[i] += step
            minus[i] -= step
            numerical.append((objective(plus, W)[0]-objective(minus, W)[0])/(2*step))
        np.testing.assert_allclose(gradient, numerical, rtol=2e-6, atol=2e-9)
        coarse = kernel.coarse_sum(x, mapping, blocks)
        np.testing.assert_allclose(coarse.sum(0), sums, rtol=1e-13, atol=1e-10)
        original, _, geometry = kernel.fine_sufficient(xy, codes, 30., 4)
        old = original["all"]
        np.testing.assert_array_equal(shape, geometry["shape"])
        assert blocks == geometry["coarse_blocks"]
        standard = kernel.spatial_estimate(old["x"], old["W"], shape, blocks)
        shared = spatial(old["x"], old["W"], shape, blocks, kernel.value_gradient)
        np.testing.assert_allclose(shared["se"], standard["se"], rtol=1e-12, atol=1e-12)
        print(setting, "point, derivatives, aggregation and shared spatial formula PASS")


if __name__ == "__main__":
    main()
