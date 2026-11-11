# Fixed robustness benchmark and completed-field comparator results

Verified nested fields 1,200/1,200. Verified paired comparisons 144/144.
The prespecified paired matrix is complete.

| Condition | Model | Nested alerts / 100 | Nested coverage / 100 | Paired nested alerts | SPINNA omnibus alerts | SPINNA improvement alerts | Verified pairs / 12 |
|---|---|---:|---:|---:|---:|---:|---:|
| 12 nm triangle | Monomer/dimer | 3 | 96 | 0 | 0 | 1 | 12 |
| 12 nm triangle | 10% trimer | 66 | 92 | 8 | 3 | 6 | 12 |
| 12 nm line | Monomer/dimer | 3 | 95 | 1 | 1 | 0 | 12 |
| 12 nm line | 10% trimer | 54 | 94 | 6 | 2 | 4 | 12 |
| 24 nm triangle | Monomer/dimer | 0 | 98 | 0 | 0 | 0 | 12 |
| 24 nm triangle | 10% trimer | 39 | 96 | 5 | 2 | 5 | 12 |
| Low retention | Monomer/dimer | 0 | 95 | 0 | 1 | 0 | 12 |
| Low retention | 10% trimer | 11 | 96 | 1 | 0 | 1 | 12 |
| Lower density and six codes | Monomer/dimer | 1 | 95 | 0 | 0 | 1 | 12 |
| Lower density and six codes | 10% trimer | 80 | 97 | 11 | 3 | 6 | 12 |
| Nonuniform intensity | Monomer/dimer | 4 | 95 | 1 | 0 | 1 | 12 |
| Nonuniform intensity | 10% trimer | 52 | 89 | 7 | 4 | 9 | 12 |

The paired sample contains 12 prespecified fields per cell. Exact binomial intervals, paired discordances, and coverage flags are in SUMMARY.json.

The comparator uses official SPINNA fits with two custom approximate fitted-null rank tests. It receives the true nuisance settings and sees positions without the permanent codes. The methods use different information and computational budgets.

The nonuniform-intensity alternative failed the prespecified coverage check. The simulations test one component model family. Experimental data with unobserved source groups require separate validation.

The original interrupted attempt remains in the provenance. Its partial restricted fit is reported separately from completed-pair costs.
