# Crosslabel and direct triple correlation on Nup96

All four scores detected all 100 constructed four-code fields at the
common reference threshold. Raw and connected direct triple correlation
each flagged 3 of 100 two-code controls. The Crosslabel point and
lower-endpoint scores each flagged 4. The comparison shows similar
decisions at this operating point. Detection reached its ceiling for
all four scores.

## Two observed catalogues

The analysis used all deposited centers, their original four codes and
one 30 nm radius, without pore assignments or further cropping.

| Source | Centers | Crosslabel estimate [95% model interval] | Raw triple correlation | Connected triple correlation | Connected contrast per anchor |
|---|---:|---:|---:|---:|---:|
| Reinhardt | 15,315 | 0.781 [0.715, 0.844] | 665.925 | 584.357 | 0.895 |
| Zaza | 13,212 | 0.553 [0.513, 0.596] | 229.828 | 187.458 | 0.620 |

The raw and connected correlation columns are dimensionless averages
over the two neighbor discs. The per-anchor conversion multiplies the
connected result by half the squared expected pooled disc count. It
puts the connected result on a descriptive count scale. Crosslabel
uses local annulus compensation and a separate pair-covariance correction.
Agreement concerns positive connected signal in both catalogues. Each
catalogue follows source picking and molecular-map processing, and
neither has a calibrated field decision in this comparison. The model
intervals describe Crosslabel. There is no paired interval for the
difference between estimators. [Reinhardt deposit](https://doi.org/10.5281/zenodo.7795826) ·
[Zaza deposit](https://doi.org/10.6084/m9.figshare.28741373)

## Paired fields made from measured localizations

The original 100 Reinhardt pairs hold positions and counts fixed while
changing code occupancy. The 99 independently generated two-code
reference fields use the same underlying template bank. Every method
has a separate upper-tail rank on the same reference fields and the
same fixed .05 cutoff. All 200 evaluation fields have supported scores.

| Score used for ranking | Two-code alerts [exact 95% interval] | Four-code alerts [exact 95% interval] | Alternative score higher within pair |
|---|---:|---:|---:|
| Raw direct triple correlation | 3/100 [0.6%, 8.5%] | 100/100 [96.4%, 100%] | 100/100 |
| Connected direct triple correlation | 3/100 [0.6%, 8.5%] | 100/100 [96.4%, 100%] | 100/100 |
| Crosslabel point estimate | 4/100 [1.1%, 9.9%] | 100/100 [96.4%, 100%] | 100/100 |
| Crosslabel lower interval endpoint | 4/100 [1.1%, 9.9%] | 100/100 [96.4%, 100%] | 100/100 |

The two triple-correlation scores flagged the same three controls.
Crosslabel's lower endpoint agreed on two, flagged two others, and
missed one alert made by both correlation scores. The exact paired
test for those three disagreements gives p = 1.00. Crosslabel's point
and lower endpoint flagged four controls each and agreed on two.
The binomial intervals are conditional on the realized reference bank.
These are template draws from one acquisition. Their code transformation
can alter lower-order summaries. The all-positive alternative result
limits conclusions about discrimination on weaker alternatives.
The earlier 2/100 false-alert result used the fixed positive-lower-limit
rule. The current 4/100 uses a common nominal reference tail.

## Comparator and replay

The independent implementation integrates the coordinate enumeration
in Yin et al. over two 30 nm discs, weights distinct code triples by
their available record counts, and applies the joint translation
correction for a rectangle. The connected version subtracts all three
integrated pair terms. The original software reports a full geometry
profile with angular edge correction and Windows MATLAB executables.
This comparison concerns the specified scalar adaptation and its shared
decision protocol. [Direct triple correlation paper](https://doi.org/10.1038/s41467-018-08048-2) ·
[Upstream software](https://github.com/yiny02/direct-Triple-Correlation-Algorithm)

The optimized counts match an independent ordered enumeration that
includes coincident records, unequal code counts and boundary cases.
All 105 package tests passed. The primary run completed in 54.5 s on
one worker. The primary Crosslabel computations were 99 reference
fits and two verification fits on the first prior pair. All 200 saved
Crosslabel scores were reused after checking the regenerated rosters,
geometry, per-code counts and witness tallies. The raw and connected
correlations share enumeration. A full second local execution
reproduced both native rows, all 99 references and all 200 evaluation
rows byte for byte in 57.0 s. The score audit independently recomputed
the ranks and paired summaries from the saved rows.
[Protocol](PROTOCOL.md) ·
[Replay receipt](REPLAY.json)

## Implication for the manuscript

Report the shared positive signal on the two measured catalogues
and the comparable decision results on the paired controls. The
result supports a comparison with a published higher-order statistic.
The contribution remains the restricted-library formulation, the
Crosslabel model interval, and the workflow for applying the check.
The earlier pair-summary result applies to its specific comparator.
The present comparison supplies no evidence of a sensitivity gain
over direct triple correlation or an experimental CD20 treatment
decision.
