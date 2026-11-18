# Shared interval for connected triple correlation

This is a retrospective ablation on all 540 saved evaluation fields from
the three-setting matched-measure follow-up. No fields, scores, cutoffs,
or settings are selected from the new interval outcomes.

Reproduce the full rectangular-translation-corrected connected triple score
per anchor. Partition its code counts, pair terms and anchor-centered
triple terms onto the same 100 nm grid used by Crosslabel. No tile-edge
neighborhoods are discarded. Treat the observed rectangular bounds as fixed
when differentiating and bootstrapping the additive contributions.

Use analytic derivatives of the complete count-normalized score, centered
fine-cell scores, the same 600 nm spatial covariance, occupied 2 micrometre
coarse blocks, 999 outer and 199 inner draws, smoothed CDF prepivot and
unchanged positive lower endpoint rule. Centering leaves Crosslabel's
homogeneous statistic unchanged, since its aggregate cell score is zero.
For the nonhomogeneous fixed-window triple statistic it removes the cell
score mean. There is no fitted variance multiplier or truth calibration.

Before collecting outcomes, test exact full-field score parity, analytic
derivatives by central differences, aggregation, and agreement of the
shared spatial formula with Crosslabel. Regenerate all 540 catalogues from
their recorded seeds and verify record count, target, Crosslabel interval
and saved triple point within documented numerical tolerances.

Report all 100 adequate and 80 alternative fields per setting, unsupported
intervals, both comparator decision rules, false alerts and paired
disagreements. The original tile-grid comparison stays in the record.
Recompute ten predetermined fields per setting after collection, including
both classes. Keep all saved sufficient arrays and nested calibration
arrays. No method repair or threshold search follows from this ablation.
