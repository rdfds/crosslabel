# Support and processing boundaries

These results are retained alongside the six-condition native all-label qualification. They prevent interpreting that pass as general processing robustness.

| Condition | View | Mask | Planned | Usable | Unsupported | Covered | Positive |
|---|---|---|---:|---:|---:|---:|---:|
| sparse_dimer | native | all | 400 | 400 | 0 | 372 | 14 |
| sparse_dimer | native | non_consecutive | 400 | 0 | 400 | 0 | 0 |
| sparse_trimer | native | all | 400 | 400 | 0 | 381 | 399 |
| sparse_trimer | native | non_consecutive | 400 | 0 | 400 | 0 | 0 |
| dense_dimer | native | all | 400 | 400 | 0 | 381 | 8 |
| dense_dimer | native | non_consecutive | 400 | 400 | 0 | 381 | 9 |
| dense_dimer | merged_10nm | all | 10 | 10 | 0 | 5 | 0 |
| dense_dimer | merged_10nm | non_consecutive | 10 | 10 | 0 | 9 | 0 |
| dense_four | native | all | 400 | 400 | 0 | 386 | 400 |
| dense_four | native | non_consecutive | 400 | 400 | 0 | 373 | 400 |
| dense_four | merged_10nm | all | 10 | 10 | 0 | 0 | 10 |
| dense_four | merged_10nm | non_consecutive | 10 | 10 | 0 | 0 | 10 |
| raw_dimer | native | all | 400 | 400 | 0 | 376 | 9 |
| raw_dimer | native | non_consecutive | 400 | 0 | 400 | 0 | 0 |
| raw_trimer | native | all | 400 | 400 | 0 | 380 | 399 |
| raw_trimer | native | non_consecutive | 400 | 0 | 400 | 0 | 0 |

The non-consecutive mask is unsupported on all400 fields in each sparse/raw condition. These are unsupported results, not successful zero calls. It supports the two dense native conditions; that secondary mask never replaces the frozen all-label primary.

The20 merged views are paired processing negatives, not additional independent observations. After10nm merging, the all-label interval covers the original target in only5/10 dense-dimer views and0/10 dense-four views. Dense-four retains10/10 positive calls despite0/10 coverage. This explicitly shows that a positive higher-order signal and the native finite-model qualification do not validate arbitrary preprocessing or identify the original biological target. The non-consecutive dense-four mask also covers0/10. No processing-robustness or biological-source claim is justified.

There were 3,240 estimated masks across the two execution segments. Shared pipeline times cover multiple estimators.
