# Held-out comparison, experimental catalogue and density screen

The design and decision rules were specified before generating the evaluation fields.
The intended library decision concerns the stated component and observation
model. The data identify where the current procedure can help and where an
experimental decision still needs a calibration study.

| Result | Observed outcome |
|---|---|
| Held-out matched-moment trimer fields | 89/100 three-code detections, 8/100 marked-pair detections |
| Same reference tail, adequate fields | 4/100 and 5/100 false alerts, respectively |
| Paired disagreements on alternatives | 82 three-code only, 1 marked-pair only, exact paired p = 1.74e-23 |
| Original Nup96 RESI center catalogue | 15,315 centers, estimate .781, descriptive model interval [.715,.844] |
| Paired controls from measured pore localizations | 100/100 detections and 2/100 false alerts with identical positions and counts within each pair |
| Density screen | 200/200 flags in the uneven-intensity condition and 3/100 in matched uniform alternatives |
| Density screen transferred to uniform lower-density six-code fields | 192/200 flags, showing the need to match reference density and sampling |

The two scores shared the same 100 adequate reference fields and nominal
tail rule. Their simulated null was specific to the matched-moment case.
The historical SPINNA comparison remains inconclusive for superiority and
is reported with its discordances and exact test. The original 89/100
coverage failure remains. The density screen flags coarse variation or a
reference mismatch. It provides no guarantee of calibration.

The original Nup96 catalogue was prepared after picking pores. The field
estimator used its x/y coordinates and codes without group IDs. The paired
fields retain measured localization geometry and repeat structure but
randomize parent centers and code occupancy. They represent one source
acquisition. They do not supply a new biological inference about CD20 or a
calibrated native-cell rejection. The published CD20 acquisition and selection
limitations are detailed in the previous report.

All 401 field calculations reproduced exactly with the packaged numerical
routines. The 200 marked-pair fits and their target counts also reproduced.
The source conversion reproduced the canonical reductions for the three
locally available original Nup96 HDF5 localization members. The full
401-call check used retained arrays. The package passed 105 tests, and
two CSV examples plus a Picasso equivalent passed from a fresh installed
wheel.

Source: [Reinhardt et al.](https://doi.org/10.1038/s41586-023-05925-9),
[localization data](https://doi.org/10.5281/zenodo.7795826).
