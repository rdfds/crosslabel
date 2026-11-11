# Triple correlation on the fixed robustness benchmark

Evaluate the same integrated direct triple correlation used on
[Nup96](../nup96/PROTOCOL.md) on all three settings of the
[original robustness study](../robustness_reference/PROTOCOL.md).
Regenerate the original 531 reference and 540 evaluation catalogues
from their pinned generator and seeds. Verify the schedule, record
counts, realized targets, and marked-pair fits against the saved rows.
Reuse the stored Crosslabel point and lower endpoint. This is a new
comparator analysis of previously evaluated fields. The benchmark
acquires no new independent experimental or simulation fields.

The primary comparison is the Crosslabel lower endpoint against
connected triple correlation at 30 nm, with the same rectangular
translation correction and connected pair subtractions as in the
Nup96 comparison. Raw triple correlation is included for context.
Also rank the Crosslabel point estimate and the connected contrast
rescaled to expected records per anchor. The rescaling varies with
each field's observed density, so this sensitivity checks whether
normalization contributes to differences in score rankings. The
point-versus-lower comparison isolates the choice of Crosslabel
score. The Crosslabel point versus direct connected comparison
changes several estimator components together and cannot isolate
the contribution of local compensation alone. Retain the original
marked-pair results on the same fields.

For each setting and each score use the same 59 reference fields at
each dimer fraction .25, .50 and .75. Take the maximum of the three
inclusive upper-tail ranks and flag at .05. If any score is
undefined in a reference group, leave that score's decisions
unresolved for the setting. Report support and exact binomial
intervals per tested composition or alternative group, with pooled
on-grid and off-grid intervals labeled descriptive. Report paired
decision disagreements against the Crosslabel lower endpoint,
including the point-versus-lower comparison. Show detection and
false alerts together in each setting. The nominal cutoff is shared
and the realized false-alert rates can differ. The grid bound applies
to represented compositions under the reference observation law.
Off-grid behavior is empirical.

Use the same independently implemented scalar triple-correlation
adaptation as in the Nup96 study. This evaluation does not rerun the
upstream MATLAB geometry profiler or compare against every possible
local-intensity adjustment. Use one worker and a 10-minute limit
for each of the primary run and one exact replay. Save every score
and verify all ranks and aggregate summaries. Do not change the
methods or thresholds in response to outcomes.
