# Connected triple correlation on the fixed robustness fields

Crosslabel's lower endpoint detected 10 more alternatives at low
retention and 15 more under varying intensity than connected direct
triple correlation at the common nominal 5% reference cutoff. It also
raised one and four additional false alerts, respectively. In the
balanced, lower-density setting the two rules detected the same 78
alternatives, with more on-grid false alerts for Crosslabel.

Every score used the original 177 reference fields per setting,
with 59 at each adequate composition. We reused the original
maximum of inclusive upper-tail ranks. The 540 evaluation fields
and all 531 references were regenerated from the original seeds.
Record counts, realized targets, fitted pair fractions and pair
scores matched the saved values on every field. Crosslabel's
original lower-endpoint decisions and marked-pair decisions were
recomputed from saved scores and matched exactly. All six scores
were supported in all three settings.
[Protocol](PROTOCOL.md) ·
[Original study](../robustness_reference/REPORT.md)

## Detection and false alerts at the fixed reference cutoff

| Setting | Score | Alternatives detected /80 | On-grid false alerts /60 | Off-grid false alerts /40 |
|---|---|---:|---:|---:|
| Balanced, lower density | Connected triple correlation | 78 | 3 | 1 |
| | Connected contrast per anchor | 78 | 4 | 1 |
| | Crosslabel point | 78 | 3 | 1 |
| | Crosslabel lower endpoint | 78 | 6 | 1 |
| | Raw triple correlation | 44 | 2 | 0 |
| | Marked pair | 0 | 2 | 0 |
| Low retention, skewed codes | Connected triple correlation | 11 | 0 | 1 |
| | Connected contrast per anchor | 12 | 0 | 1 |
| | Crosslabel point | 14 | 0 | 1 |
| | Crosslabel lower endpoint | 21 | 1 | 1 |
| | Raw triple correlation | 1 | 0 | 0 |
| | Marked pair | 3 | 3 | 1 |
| Varying intensity, skewed codes | Connected triple correlation | 27 | 0 | 0 |
| | Connected contrast per anchor | 29 | 0 | 0 |
| | Crosslabel point | 34 | 0 | 1 |
| | Crosslabel lower endpoint | 42 | 1 | 3 |
| | Raw triple correlation | 1 | 0 | 0 |
| | Marked pair | 0 | 0 | 0 |

The two primary scores detected the same 78 fields in the balanced
setting. At low retention and under varying intensity, the additional
10 and 15 detections were all Crosslabel-only. The two-sided exact
paired p-values on alternatives are .001953 and .0000610,
respectively, conditional on the reference banks. These tests compare
the decisions at the fixed nominal cutoff. The realized false-alert
rates differ, so these comparisons do not establish detection gains
at an equal realized false-alert rate. The per-composition exact
intervals, paired disagreements, support and all other comparisons
are in [RATES.csv](RATES.csv). Pooled on-grid and off-grid intervals
are descriptive because the fields come from different compositions.

## What the components add

Rescaling the connected triple correlation to records per anchor
increased detections from 11 to 12 at low retention and from 27 to
29 under varying intensity. Crosslabel's point estimate detected
14 and 34, respectively. The point comparison changes several
elements together, including local compensation, the pair correction
and normalization. These results do not identify which one caused
the differences.

Ranking the Crosslabel lower endpoint added seven and eight
detections relative to its point estimate in those two settings.
It also added one on-grid alert at low retention and one on-grid
and two off-grid alerts under varying intensity. In the balanced
setting it added three on-grid alerts with unchanged detections.
The point and lower endpoint use the same estimator, so this
comparison isolates the choice of score. Interval coverage is a
separate question. The alternative coverage counts from the
original study remain 74, 77 and 75 of 80 in setting order.

The biological comparison shows positive connected signals in
both public Nup96 catalogues and 100/100 detections for each rule
on the paired templates. Together, the results support a restricted
library workflow with a model-based uncertainty score that changes
decisions under specific simulated observation laws. They leave
experimental calibration and matched-error superiority open.
The comparator is the same 30 nm integrated adaptation of direct
triple correlation used on Nup96, with rectangular translation
correction and all three pair terms subtracted.
[Nup96 comparison](../nup96/REPORT.md) ·
[Direct triple correlation paper](https://doi.org/10.1038/s41467-018-08048-2)

The complete replay reproduced both JSONL files byte for byte.
All 3,240 saved evaluation ranks and all 144 rate rows also
recompute from the saved records. The replay regenerates the same
coordinates and does not add independent fields.
[Replay receipt](REPLAY.json) · [Full summary](SUMMARY.json)
