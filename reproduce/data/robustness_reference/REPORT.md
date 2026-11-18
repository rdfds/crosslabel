# Robustness and comparator simulation study

The fixed experiment evaluated the existing Crosslabel lower limit and marked-pair
residual on the same simulated catalogues. Both scores used the same prespecified
reference grid and upper-tail decision rule within each observation setting.
The three new settings combine changes to retained density, code probabilities,
retention and parent intensity. The protocol was specified before generating
the reference and evaluation fields.

| Setting | Crosslabel detections (95% CI) | Marked-pair detections (95% CI) | Crosslabel only | Pair only | Model interval coverage |
|---|---:|---:|---:|---:|---:|
| Balanced codes, lower density | 78/80 (91.3–99.7%) | 0/80 (0.0–4.5%) | 78 | 0 | 74/80 |
| Low retention, skewed codes | 21/80 (17.0–37.3%) | 3/80 (0.8–10.6%) | 19 | 1 | 77/80 |
| Varying intensity, skewed codes | 42/80 (41.0–63.8%) | 0/80 (0.0–4.5%) | 42 | 0 | 75/80 |

| Setting | On-grid Crosslabel false alerts | On-grid pair false alerts | Off-grid Crosslabel false alerts | Off-grid pair false alerts |
|---|---:|---:|---:|---:|
| Balanced codes, lower density | 6/60 (3.8–20.5%) | 2/60 (0.4–11.5%) | 1/40 (0.1–13.2%) | 0/40 (0.0–8.8%) |
| Low retention, skewed codes | 1/60 (0.0–8.9%) | 3/60 (1.0–13.9%) | 1/40 (0.1–13.2%) | 1/40 (0.1–13.2%) |
| Varying intensity, skewed codes | 1/60 (0.0–8.9%) | 0/60 (0.0–6.0%) | 3/40 (1.6–20.4%) | 0/40 (0.0–8.8%) |

All 540 evaluation calls had an interval and a reference rank. Detection
limits are exact binomial intervals conditional on the realized reference
banks. The false-alert table pools fixed composition counts within each
stratum and uses descriptive binomial limits. The true rates can vary by
composition. The JSON and CSV retain exact limits for each separate
composition and exact paired tests within each setting.

The balanced lower-density setting shows a detection gain alongside more
false alerts at the evaluated threshold. Low retention still limits detection
to 21 of 80 alternatives despite improvement over the tested pair score.
Each setting varies several nuisance inputs together. The reference procedure
assumes the specified observation law. Its finite-grid bound relies on
exchangeability of the reference and test draws. Experimental decisions require
independent characterization and calibration of the acquisition.

Rank calibration addresses the score decision under the stipulated simulations.
The unmodified interval coverage remains a separate quantity. Earlier undercoverage
under uneven intensity, the density-screen transfer failure, and the inconclusive
SPINNA comparison remain part of the evidential record. The simulation settings
and reference thresholds were specified separately from the experimental results.

Collection completed 1,071 reference and evaluation calls and 384 pair-only
training fields in 352.4 seconds. The saved per-field
scores, intervals, targets and training histograms accompany the summary.

A clean archive extraction and installed Crosslabel 0.4.0 wheel replayed all
1,071 field calls and 384 pair-training histograms exactly. The replay
used the same numerical dependency versions and took 405.8 seconds.
Both installed CLI examples matched their expected reports. The archived
earlier 200 held-out marked-pair fits also reproduced. The archive loader
verified all eight retained Nup96 arrays and reconstructed the recorded
roster of 15,315 centers and 1,206 templates. The package suite passed
105 tests.
