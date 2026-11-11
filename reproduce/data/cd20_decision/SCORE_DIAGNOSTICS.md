# Descriptive comparison of saved scores

The prespecified pair rank was conservative in the completed study.
This post hoc diagnostic compares the ordering given by the unchanged
pair residual and the lower endpoint of the nested interval.
Both scores use the same fields. The table describes the empirical ROC
at 5 false alerts among 100 simulated adequate fields and the mean
sensitivity over the observed 0–5% false-alert range. Ties use
linear interpolation. This is an in-sample descriptive ROC comparison.
It supplies neither a new calibrated threshold nor confirmation.

| Family | Nested sensitivity at 5% | Pair sensitivity at 5% | Nested mean sensitivity 0–5% | Pair mean sensitivity 0–5% |
|---|---:|---:|---:|---:|
| reference | 80% | 4% | 55.2% | 1.4% |
| linear | 65% | 5% | 41.0% | 2.4% |
| wide | 53% | 2% | 44.8% | 2.0% |
| low_retention | 35% | 4% | 17.6% | 2.4% |
| low_density_six_codes | 94% | 0% | 86.0% | 0.0% |
| known_intensity | 59% | 1% | 44.4% | 0.4% |
| matched_moments | 96% | 4% | 92.4% | 3.2% |
| membrane_variation | 51% | 2% | 37.6% | 1.2% |
