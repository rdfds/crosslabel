# Field area at fixed retained density

Use the same homogeneous triangular simulator as the retention study.
Retained subunit density stays 70 per square micrometre. Compare nominal
parent-window areas 100, 200 and 400 square micrometres, at retention
0.25, 0.50 and 0.75. Expected observed anchors are 7,000, 14,000 and 28,000.
The observation error, four balanced codes, fractions (0.5,0.5,0) or
(0.7,0.2,0.1), 30 nm radius and Crosslabel interval are fixed.

Reuse all 80 adequate and 80 alternative 100-square-micrometre fields at
each retention from the completed fixed-count study. This is 480 reused
fields. Generate exactly 960 additional independent fields for the two
larger areas, with 80 per class and setting. New SeedSequence entropy
is [19112026, area_index, retention_index, alternative, replicate],
where area_index is 1 or 2, retention_index is 0,1,2, and replicate 1..80.
Bootstrap seed is 19000000 + 100000*area_index + 10000*retention_index
+ 1000*alternative + replicate. No replacements or sample expansion.

Scale the Poisson mean by area and sample centers over the stated larger
square. Reproduce every baseline catalogue and target from the prior
schedule before using the generalized generator. Keep all sufficient raw
catalogues and new observations. Bootstrap blocks retain the same physical
2 micrometre side; the covariance radius remains 600 nm.

Report all nine settings and both classes, realized counts, detection,
false alerts, coverage, unsupported intervals and exact binomial limits.
The analysis is a supplementary sensitivity study, with previously seen
baseline outcomes. It is separate from the six-condition interval
qualification and does not tune any method or threshold.

One 45-minute collector with an exclusive lock and 2 GiB reserve. An
interruption preserves evidence and does not trigger an automatic retry.
Regenerate and fully replay the first adequate and alternative field
in every new cell after collection.
