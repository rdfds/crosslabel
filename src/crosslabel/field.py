"""Field-level connected three-code contrast and approximate interval."""
from ._field_kernel import (
    coarse_sum, fine_sufficient, nested_estimate, spatial_estimate, validate_catalogue,
)
from .field_screen import density_score, reference_tail


def field_check(xy, codes, *, code_count, radius_nm, seed, reference_scores=None):
    """Check one complete 2D catalogue with 4 or 6 permanent codes.

    Coordinates are in nm and code IDs are integers starting at zero. Use 10
    or 30 nm radii as in the validated workflows. The optional density
    reference must come from independently justified comparable fields.
    Results describe the specified observation and component model. Using
    the interval for an experimental decision requires separate calibration.
    """
    xy, codes, support = validate_catalogue(xy, codes, code_count, radius_nm, "nm")
    if isinstance(seed, bool) or int(seed) != seed or seed < 0:
        raise ValueError("A nonnegative integer random seed is required")
    seed = int(seed)
    density = density_score(xy)
    density["reference_p"] = (
        None if reference_scores is None else reference_tail(density["score"], reference_scores)
    )
    density["flagged"] = None if reference_scores is None else density["reference_p"] <= .05
    masks, mapping, geometry = fine_sufficient(xy, codes, radius_nm, code_count)
    all_codes = masks["all"]
    blocks = geometry["coarse_blocks"]
    if blocks != support["coarse_blocks"]:
        raise ValueError("Input support differs from evaluated geometry")
    coarse = coarse_sum(all_codes["x"], mapping, blocks)
    spatial = spatial_estimate(all_codes["x"], all_codes["W"], geometry["shape"], blocks)
    result, _ = nested_estimate(coarse, all_codes["W"], spatial, seed)
    interval = result.get("spatial_bootstrap_t_ci")
    return {
        "status": result["status"],
        "reason": result.get("reason", spatial.get("reason")),
        "estimate": result.get("value"),
        "interval_95_model": interval,
        "positive_lower_limit": None if interval is None else interval[0] > 0,
        "density_screen": density,
        "support": support,
        "window_nm": geometry["window"],
        "seed": seed,
        "radius_nm": radius_nm,
        "interpretation": "model_based_contrast_requires_applicability_review",
    }
