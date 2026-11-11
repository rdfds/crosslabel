# Field-only decisions and CD20 catalogue contrasts

We applied both rules to the same 540 previously saved evaluation fields.
Crosslabel flags a positive lower endpoint of its nominal 95% model
interval. The direct triple-correlation comparator divides the observed
rectangle into 25 equal tiles, computes the connected score separately in
each, and flags a positive lower endpoint of a nominal 95% Student interval
for the mean tile score. Both rules use only the evaluation field. The
published score by itself has no decision threshold in our implementation.
The tile rule is one approximate way to supply uncertainty without a
simulated reference bank. It targets the tile mean and treats tiles as
approximately independent and comparable, even though neighboring tiles
may be dependent and spatial intensity may vary.

| Setting | Crosslabel false alerts /100 | Tiled triple correlation false alerts /100 | Crosslabel detections /80 | Tiled triple correlation detections /80 |
|---|---:|---:|---:|---:|
| Balanced, lower density | 5 | 0 | 78 | 76 |
| Low retention, skewed codes | 2 | 0 | 22 | 8 |
| Varying intensity, skewed codes | 4 | 4 | 46 | 31 |
| Uniform-intensity control matched to varying setting | 3 | 1 | 54 | 43 |

Both rules returned a decision on every field. In the varying-intensity
setting, 17 alternatives were detected only by Crosslabel, two only by
tiled triple correlation, and 29 by both. Two of the four false alerts
were shared. The exact 95% binomial interval for 4/100 is 1.1% to 9.9%.
The observed false-alert counts are too imprecise to establish a 5% upper
bound. The low-retention gain accompanies two extra false alerts. The
comparison concerns these two decision procedures, and the simulated
reference bank results answer a different operating-point question.
[Matched false-alert comparison](../robustness/MATCHED_FALSE_ALERTS.md)

### Tile-size sensitivity

All 540 original 5 by 5 scores reproduced exactly. Adjacent grids also
returned decisions on every field. Each entry is detections per 80
alternatives, followed by false alerts per 100 adequate fields.

| Setting | Crosslabel | 4 by 4 | 5 by 5 | 6 by 6 |
|---|---:|---:|---:|---:|
| Balanced, lower density | 78, 5 | 75, 1 | 76, 0 | 72, 0 |
| Low retention, skewed codes | 22, 2 | 10, 0 | 8, 0 | 9, 0 |
| Varying intensity, skewed codes | 46, 4 | 47, 3 | 31, 4 | 30, 0 |

The coarser tile grid detected one more varying-intensity alternative than
Crosslabel and produced one fewer observed false alert. The two methods can
both make field-only decisions, and the observed comparison depends on tile
size. [All tile scores](TILE_SENSITIVITY.json)

## What changed when we removed the intensity pattern

We regenerated the 180 varying-intensity fields with uniform probabilities
over the same parent grid. The generator used the same random streams.
The record counts, marks and realized target agreed in every pair.
Replacing the intensity pattern raised Crosslabel detection from 46 to
54 of 80 and tiled triple-correlation detection from 31 to 43 of 80.
Crosslabel model-interval coverage increased from 92/100 to 95/100 on
adequate fields and from 75/80 to 77/80 on alternatives. These are modest
coverage shifts for one simulated pattern.

The Crosslabel point estimate's mean signed error on alternatives was
0.000435 under varying intensity and 0.000468 under uniform intensity.
Median total interval width was 0.0613 and 0.0571, respectively. The
nearly unchanged mean error and wider intervals are consistent with a
precision penalty, although changes in the field-by-field estimates and
intervals can both affect calls. In the separate matched-alert benchmark
on the varying-intensity fields, Crosslabel's point score detected 62/80,
its lower endpoint detected 49/80, and full-field connected triple
correlation detected 61/80 at five observed false alerts each. This
identifies the uncertainty-adjusted ranking as the source of the observed
performance drop relative to Crosslabel's own point ranking. It does not
establish which component of the interval construction drives the change.

## Descriptive obinutuzumab CD20 application

We scored all retained centres in both verified round-resolved
obinutuzumab acquisitions at 30 nm using lateral coordinates and the
published round assignments. The G5M and gradient-ascent (GA) columns
are two reconstructions of the same acquisitions.

| Acquisition | G5M centres | G5M contrast | GA centres | GA contrast |
|---|---:|---:|---:|---:|
| cell1_2023_Pos3 | 17,679 | 0.3314 | 15,075 | 0.2295 |
| cell4_3D, lateral projection | 27,167 | 0.2566 | 21,807 | 0.2057 |

The connected triple-correlation contrast per anchor was also larger for G5M
on both acquisitions (0.5048 vs 0.3566 and 0.4089 vs 0.3210).
We drew 32 uniform subsamples of G5M without replacement on each cell to
match the per-round GA counts. The median Crosslabel contrast was 0.2414
(range 0.2046–0.2838) on the first and 0.1627 (0.1380–0.1802) on the second.
The subsampled contrasts exceeded full GA on 22 and 0 draws, respectively.
The second cell's ordering reversed with count matching. The draws share
their source catalogue and do not measure between-cell variation. The
effective footprint remains unverified, and the recovered centers differ
between reconstructions. These results describe the selected catalogues.
The G5M study reports higher oligomerization
for the inverted T-cell engager relative to obinutuzumab and the
classical engager. Our verified local per-round inputs cover two
obinutuzumab acquisitions. The retained Pachmayr engager catalogues are
pooled and lack verified per-record round assignments. The available
comparison therefore concerns reconstruction on two cells from one
treatment. [G5M study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12976027/)
· [G5M deposit](https://zenodo.org/records/18429674)
· [Pooled data audit](../source_catalogues/REPORT.md)

## Reproduction

[Protocol](PROTOCOL.md)
· [Analysis](../../field_only.py)
· [Field results](FIELD_RESULTS.json) · [CD20 results](CD20.json)
· [CD20 sensitivity](CD20_SENSITIVITY.json)

Run from the repository root with the validation Python environment.
The script verifies the regenerated original fields and the CD20 input
checksums. Use fresh output paths because the script opens them exclusively.

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B \
  reproduce/field_only.py /tmp/crosslabel-field-results.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B \
  reproduce/tile_sensitivity.py /tmp/crosslabel-tiles.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B \
  reproduce/field_only.py --cd20-sensitivity --g5m-dir /path/to/extracted-g5m /tmp/crosslabel-cd20.json
```
