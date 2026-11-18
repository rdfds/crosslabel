# Retention, observed anchor count and execution cost

This fixed supplementary study was specified before its outcomes. It uses
the released Crosslabel field interval, without a reference bank or changes
to its estimator, bootstrap draws, covariance radius or decision threshold.
It supplements the earlier independently seeded qualification and does not
replace or pool that study.

Use the existing compound-Poisson generator in
generator.py, with homogeneous independent parent
centers in a 10 by 10 micrometre window, 12 nm triangular geometry, the
original positional noise, four equiprobable permanent source codes,
one record per retained source, and a 30 nm anchor radius.
The retained subunit density is the expected anchor count divided by
100 square micrometres. Thus retention changes while expected observed
count stays fixed within each curve. This holds observed density fixed
by increasing the unobserved parent intensity as retention decreases.
The resulting curves apply to this simulator family and footprint.

Cross retention probabilities 0.15, 0.25, 0.50, 0.75 and 1.00 with
expected observed record counts 3500, 7000 and 14000. At each of the
15 settings generate 80 adequate and 80 alternative independent fields.
Adequate subunit fractions are (0.5, 0.5, 0); alternative fractions are
(0.7, 0.2, 0.1). There are 2400 fields. Order by replicate, retention,
count and class. RNG SeedSequence entropy is
[18112026, retention_index, count_index, alternative, replicate];
bootstrap seeds are 18000000 + 100000*retention_index +
10000*count_index + 1000*alternative + replicate (replicates 1..80).
All fields are new supplementary sensitivity data.

Flag only if the unchanged 95% interval lower endpoint exceeds zero.
Retain every unsupported case and error. Report detections and false alerts
over all 80 planned fields per cell, supported counts, target coverage,
exact binomial limits for alert proportions, realized record counts,
target and estimate distributions. No favorable stopping, parameter
selection, new repair or replacement seed follows from these outcomes.
An unsupported case is separately reported, not classified as a successful
non-alert. Only a complete supported cell receives an unqualified rate.

Save each complete catalogue, its raw hash and original seed for exact
regeneration. Regenerate the first adequate and alternative in every cell
and compare every coordinate, code, target, endpoint and decision before
reporting. This is additional execution, not independent evidence.

Measure elapsed wall time and process CPU time for the public field_check
call on every field, excluding simulation, file I/O, and process startup.
The public call includes input checks, score computation, covariance and
all 999 by 199 bootstrap draws. Report median and interquartile ranges by
expected count with hardware, OS, Python and dependency versions and
single-thread settings. Do not infer costs for individual shared-pipeline
estimators or extrapolate a SPINNA/runtime superiority comparison.

One collector, exclusive lock and original source hash checks. A 40-minute
execution bound and 2 GiB disk reserve apply. An interruption preserves
all files and stops the study without an automatic retry.
