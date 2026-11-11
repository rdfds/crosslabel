# Marked-pair comparator and CD20 application decision

**The ungrouped CD20 application remains unqualified.** All six source fields
have recorded round identities. The effective common membrane footprint and
the observation law after roundwise selection and clustering remain unresolved.
Neither a fitted fraction nor a rejection was computed on the CD20 inputs.

Collector status is **complete**. Pair scores and decisions were replayed
from saved summaries for 2400 rows. For the fresh fields we checked the saved
interval decisions and coverage against their stored targets. This verification
did not repeat nested resampling or regenerate the fresh fields.

## Fixed-threshold comparisons

Null alerts are false alerts under the stated simulated library. Alternative
alerts are detections. Discordances count only fields with a nested interval.
The first six families reuse the existing fields. The last two use independent
seeds. A field without an alert leaves the library unresolved.

| Family | Class | Saved / planned | Nested intervals | Nested alerts | Marked-pair alerts | Nested only | Pair only | Covered |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| reference | null | 100/100 | 100 | 3 | 0 | 3 | 0 | 96/100 |
| reference | alternative | 100/100 | 100 | 66 | 0 | 66 | 0 | 92/100 |
| linear | null | 100/100 | 100 | 3 | 0 | 3 | 0 | 95/100 |
| linear | alternative | 100/100 | 100 | 54 | 0 | 54 | 0 | 94/100 |
| wide | null | 100/100 | 100 | 0 | 0 | 0 | 0 | 98/100 |
| wide | alternative | 100/100 | 100 | 39 | 0 | 39 | 0 | 96/100 |
| low_retention | null | 100/100 | 100 | 0 | 0 | 0 | 0 | 95/100 |
| low_retention | alternative | 100/100 | 100 | 11 | 0 | 11 | 0 | 96/100 |
| low_density_six_codes | null | 100/100 | 100 | 1 | 0 | 1 | 0 | 95/100 |
| low_density_six_codes | alternative | 100/100 | 100 | 80 | 0 | 80 | 0 | 97/100 |
| known_intensity | null | 100/100 | 100 | 4 | 0 | 4 | 0 | 95/100 |
| known_intensity | alternative | 100/100 | 100 | 52 | 0 | 52 | 0 | 89/100 |
| matched_moments | null | 100/100 | 100 | 2 | 0 | 2 | 0 | 96/100 |
| matched_moments | alternative | 100/100 | 100 | 91 | 0 | 91 | 0 | 95/100 |
| membrane_variation | null | 100/100 | 100 | 2 | 0 | 2 | 0 | 94/100 |
| membrane_variation | alternative | 100/100 | 100 | 31 | 0 | 31 | 0 | 95/100 |

## Independent off-grid adequate controls

The rank test is calibrated for each grid composition under the specified
generator. The checks below measure off-grid performance at two compositions.
The maximum rank over the grid can be conservative. The two procedures use
different scores and thresholds. Sensitivity comparisons at an equal realized
false-alert rate would require separate operating-point calibration.

| Family | q = .25 | q = .75 |
|---|---:|---:|
| reference | 0/50 | 0/50 |
| linear | 0/50 | 1/50 |
| wide | 0/50 | 0/50 |
| low_retention | 0/50 | 0/50 |
| low_density_six_codes | 0/50 | 0/50 |
| known_intensity | 0/50 | 0/50 |
| matched_moments | 0/50 | 0/50 |
| membrane_variation | 0/50 | 0/50 |

## Six experimental fields

| Field | Treatment | Localization records | Map rows | Max/min round map count | Decision |
|---|---|---:|---:|---:|---|
| fig4untreated | untreated | 4,912,617 | 21,515 | 4.27 | Abstain before fitting |
| edf10rep1 | untreated | 3,749,188 | 69,091 | 1.81 | Abstain before fitting |
| edf10rep2 | untreated | 2,065,510 | 19,755 | 1.64 | Abstain before fitting |
| rtx_rep1 | rituximab | 3,098,768 | 20,684 | 1.40 | Abstain before fitting |
| rtx_rep2 | rituximab | 3,200,839 | 22,898 | 1.93 | Abstain before fitting |
| rtx_rep3 | rituximab | 2,512,460 | 33,998 | 1.31 | Abstain before fitting |

All six fields have four identified acquisition rounds. The source paper
documents cellular fiducial alignment, round-specific clustering and temporal
filtering. Camera ROI and coordinate bounding boxes do not specify the effective
cell footprint or how a source contributes to each round after processing.
The source already reports higher-order CD20 organization after treatment.
[RESI methods and CD20 results](https://doi.org/10.1038/s41586-023-05925-9)

The closest remaining experimental input is a field-linked mask and acquisition
or processing record that constrains spatial detection and cross-round source
assignment. Use independent controls or justified error ranges to calibrate
both rules before reading the six experimental decisions.

## Execution and reproduction

Collector wall span was 178.976 seconds, covering training,
calibration, catalogue reads, both pipelines on the fresh fields and saving.
It is a study-level measurement. Pair score replay adds verification
executions, while the 1,200 saved catalogues retain their existing nested
and SPINNA outcomes. Detailed counts are in the summary.

The reported fraction projection is a nuisance fit within a fixed simulated
library. Its finite-grid rank test has no continuous-composition guarantee.
The fresh controls use one retained record per source and stylized spatial
intensity. They do not reproduce the actual CD20 selection and clustering.

## Interpretation for the paper

The independent matched-measure control is the clearest increment over
the tested pair summary. Its two populations share first- and second-order
marked factorial measures. The fixed three-code rule detected 91/100
alternatives with 2/100 false alerts. The pair rule made no calls in either
class. A separate descriptive analysis of the unchanged scores gave 96/100
and 4/100 detections at the same observed 5/100 false-alert rate. The
comparison uses the observed fields to read the ROC and does not establish
a prospective operating threshold. See [the score comparison](SCORE_DIAGNOSTICS.md).

In the new two-dimensional intensity control, coverage was 94/100 for
adequate fields and 95/100 for alternatives. Both classes met the original
compatibility criterion, but detection reached 31/100. The prior
one-dimensional intensity alternative still failed at 89/100 coverage.
Neither intensity pattern qualifies inference for the CD20 acquisitions.

The six original combined-map HDF5 member hashes match the earlier
archive receipts. Their headers provide one localization table per member
and no root or table attributes. The recorded camera ROI differs from the
effective cellular footprint. The source paper mentions a final round with
all four imagers, which is a possible source of a coarse membrane mask
if the registered field-level record and coverage are available. This audit
did not inspect unretained archive members or infer a mask from the fitted
catalogue. [RESI source archive](https://doi.org/10.5281/zenodo.7795826).

The simulation results improve the evidence for a distinct library check.
The biological decision and a reviewer-accessible research archive remain
open for submission.

A separate deterministic replay regenerated all 400 fresh fields and
reproduced their marked-pair summaries, target counts, and saved nested
endpoints exactly. It took 151.674 seconds and adds no
independent fields. Five focused comparator tests and all 103 package
tests passed.
