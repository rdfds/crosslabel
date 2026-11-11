# Detection at matched observed false-alert counts

Across the saved fields, the Crosslabel lower endpoint tied connected triple
correlation in the balanced setting, detected one more low-retention
alternative and detected twelve fewer varying-intensity alternatives at five
observed false alerts per 100 adequate fields. Crosslabel's point score and
direct connected triple correlation had similar detection in the
varying-intensity setting at that operating point.

## Five observed false alerts per 100 adequate fields

| Setting | Connected triple correlation | Connected per anchor | Crosslabel point | Crosslabel lower endpoint |
|---|---:|---:|---:|---:|
| Balanced, lower density | 78/80 | 79/80 | 79/80 | 78/80 |
| Low retention, skewed codes | 26/80 | 23/80 | 18/80 | 27/80 |
| Varying intensity, skewed codes | 61/80 | 57/80 | 62/80 | 49/80 |

Each setting uses 60 on-grid and 40 off-grid adequate fields, with 20
fields at each of five dimer fractions. At five pooled false alerts, direct
connected triple correlation and the Crosslabel lower endpoint each flagged
four on-grid fields and one off-grid field in the balanced setting, and
three on-grid and two off-grid in the low-retention setting. Under varying
intensity the splits were two on-grid and three off-grid for direct connected
triple correlation, and one on-grid and four off-grid for the lower endpoint.
The Crosslabel point had two on-grid and three off-grid alerts in that
setting.

## Method and interpretation

For each score and setting, sweep upper thresholds on the saved field
scores and take the least restrictive threshold within each observed
false-alert budget. Include ties together. The sweep covers budgets of
zero through ten false alerts per 100 adequate fields. Each score attained
exactly five at the operating point shown above. Using the saved reference
ranks instead gave 78 vs. 73, 25 vs. 26, and 60 vs. 42 detections for direct
connected triple correlation versus the Crosslabel lower endpoint, with
actual false-alert counts of 5 vs. 0, 5 vs. 5, and 5 vs. 4. Discrete
reference ranks can skip the desired budget.

The [operating-point table](OPERATING_POINTS.csv) contains both orderings
and all 264 rows, including the realized on-grid and off-grid counts.

These operating points are descriptive. The same adequate fields determine
the thresholds and supply the observed false-alert counts. The pooled rate
refers to the evaluated mixture of five compositions. The original
prespecified nominal 5% comparisons remain in the [primary report](REPORT.md).
Use both results when interpreting the apparent detection gains. Across
these tests, which score detects more depends on the observation setting and
operating point.

A prospective comparison would set separate score thresholds on independent
adequate controls at a common target false-alert rate, then assess both
detection and false alerts on a separate paired evaluation set. For an
experimental library decision, the same design also needs independently
justified field support and observation assumptions. The [CD20 applicability
audit](../cd20_decision/REPORT.md) identifies the unresolved common
footprint and roundwise observation process in the existing fields.

Reproduce from the repository root:

```sh
python3 -B reproduce/match_false_alerts.py reproduce/data/robustness
```
