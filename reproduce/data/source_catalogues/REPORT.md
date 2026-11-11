# Public RESI data and CD20 round provenance

## Pachmayr CD20 deposit

Both archives were downloaded from [Zenodo 15552356](https://doi.org/10.5281/zenodo.15552356).
Their MD5 values match the deposit, and all 108 archive entries tested by
`lsar -t` passed. [Pachmayr et al.](https://doi.org/10.1038/s41467-025-61893-w)
describe four stochastic CD20 codes and a separate set for the antibodies.
The fitted proportions in a representative OBZ-treated cell are estimates
from the authors' model. They do not provide independent oligomer truth.

| Archive | Bytes | Zenodo MD5 | SHA-256 |
|---|---:|---|---|
| Data.rar | 14,334,971 | `31f99eced18b52a67678d7579dc7452a` | `dc35326edfce9f050441cbb4f9c14d53f5b43cd3647b12dd63a61edefae8278e` |
| Code.rar | 11,413,601 | `2d6bc40937aad8e0a2a0e10275217654` | `f93d9d638d9259b37b4686a5b08c0891efe165679ef6c707ff6304eb25ad49bf` |

Each CD20 table in Data.rar has one molecule-center dataset. The associated
sidecar identifies four source filenames under Picasso Render Combine. The
HDF5 schema contains neither a round column nor a code attribute. The
available identifiers are group and group_input. The number of distinct
group_input values varies from one to thirteen across tables. The frame
values are sorted and cover overlapping within-movie frame ranges. The
historical Picasso Combine implementation stacks and sorts the input tables
by frame and records their paths in the sidecar without adding channel IDs.
The listed filenames, group identifiers and row order therefore do not
provide a verified mapping from deposited rows to the four codes.
[Picasso Render source](https://github.com/jungmannlab/picasso/blob/v0.7.5/picasso/gui/render.py#L11335-L11349)

| Data.rar CD20 table | Centers | Distinct group_input values |
|---|---:|---:|
| 2H7 | 38,279 | 1 |
| c-TCE | 6,518 | 12 |
| H299 | 5,433 | 5 |
| i-TCE | 12,687 | 13 |
| OBZ_2D | 13,549 | 7 |
| OBZ_3D | 24,475 | 7 |
| OFA | 27,926 | 1 |
| RTX_2D | 11,163 | 8 |
| RTX_3D | 19,699 | 10 |

The separate Code.rar contains pooled center tables, downstream DBSCAN
tables and simulation files. The public scripts define the binary
`protein` column as CD20 versus antibody. Neither archive supplies a
verified per-row CD20 round map for an OBZ versus RTX comparison. No
author-provided scripts were executed.

## G5M cross-check

The existing G5M member index identifies per-round OBZ inputs. The four
metadata members inspected in five bounded requests are recorded in
[the receipt](g5m_receipt.json). The original R1 movie basenames agree
exactly between Pachmayr OBZ_2D and G5M cell1_2023_Pos3 and between
Pachmayr OBZ_3D and G5M cell4_3D. Treat each link as reuse of the
same source acquisition. The checked G5M sidecars record circular
pick diameters and a null camera ROI entry. They supply no coordinates
for a common effective footprint or independent bound on selection
and clustering. The previous [acquisition audit](../acquisition_audit/REPORT.md)
already checked the two G5M cells and their per-round files.
[G5M paper](https://doi.org/10.1038/s41467-026-70198-5) ·
[G5M deposit](https://doi.org/10.5281/zenodo.18429674)

## Second-laboratory demonstration

The public Zaza Nup96 Fig. 3 table contains 13,212 centers, with original
round counts 2,638, 3,538, 3,364 and 3,672. The extracted HDF5 and YAML
both passed ZIP member CRC. The sidecar records 108 nm per pixel and a
polygon pick. Crosslabel ran on all deposited centers at the prespecified
30 nm radius without reading any group or cluster labels. Its estimate
was 0.553 and the 95% model interval was [0.513, 0.596]. This is one
catalogue from a different microscope and laboratory than the Reinhardt
Nup96 example. The interval describes the fitted model. Its coverage
and library decision have not been calibrated for this field. There was
no CD20 model fit, treatment decision, or new comparator in this audit.
[Zaza et al.](https://doi.org/10.1038/s41467-025-60263-w) ·
[Figshare data](https://doi.org/10.6084/m9.figshare.28741373) ·
[Exact result](ZAZA_DESCRIPTIVE.json)
