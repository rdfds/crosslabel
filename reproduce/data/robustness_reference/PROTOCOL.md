# Robustness of the comparison at a common reference tail

This follow-up asks whether the advantage over marked pairs transfers across
three fixed observation settings. It extends the single-composition,
single-setting result using the existing independent Poisson component
generator and unchanged Crosslabel interval. The paired scores are the
Crosslabel lower endpoint and the fitted marked-pair residual.

All settings use four codes, a 10 micrometre square for component centers,
7 nm independent Gaussian member displacements per coordinate, 5.025 nm
independent localization error per coordinate and a 30 nm analysis radius.
Use the exact existing localization error from the generator.

| Setting | Expected retained records per square micrometre | Retention | Code probabilities | Center intensity |
|---|---:|---:|---|---|
| Balanced, lower density | 38 | .60 | .25, .25, .25, .25 | Uniform |
| Low retention, skewed codes | 55 | .27 | .10, .20, .30, .40 | Uniform |
| Spatially varying, skewed codes | 65 | .47 | .15, .26, .48, .11 | Original fixed two-dimensional mask |

Train the marked-pair fit on 64 monomer-only and 64 dimer-only fields per
setting. The pair statistic uses all code counts and pair distance bins, with
the original covariance shrinkage and projection. Training, calibration
and evaluation use disjoint fixed seed streams.

For each setting and each dimer subunit fraction q in {.25,.50,.75}, draw
59 adequate reference fields. For each score take its inclusive Monte
Carlo upper-tail rank within each composition and then the maximum of
the three ranks. Flag at .05. At a represented composition, the maximum
is at least the valid rank from the correct reference law. This yields a
finite-grid level bound conditional on the pair training and exact
simulation law, with the usual Monte Carlo exchangeability over the
reference and test draws. It does not imply a conditional level bound
given the realized reference bank, validity between grid points, or
validity under mismatched nuisance inputs.

Evaluate 20 independent adequate fields at each represented composition
and 20 each at the intervening compositions .375 and .625 per setting.
Evaluate 80 independent alternative fields per setting, with monomer,
dimer and trimer subunit fractions .65, .20 and .15. The alternative
matches the first- and second-order marked factorial measures of the
q=.50 restricted library because total subunit intensity is held fixed,
members have identical independent displacement and observation laws,
and .20 + 2(.15) = .50. The third factorial measure differs. This
matching concerns population measures, which allows finite-field pair
scores to differ through higher-order variation.

Report every setting and composition, exact binomial intervals for
detection and false alerts, paired discordance, unsupported calls and
unmodified interval coverage against the simulated realized target.
No averaging across settings determines an overall pass. A positive
comparison requires attention to false-alert tradeoffs. Keep the earlier
single-composition result, comparisons against SPINNA, and observed
coverage failures. No claim of experimental calibration follows.
