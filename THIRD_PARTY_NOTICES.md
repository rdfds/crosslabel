# Third-party data and software

Crosslabel software is distributed under the MIT license in `LICENSE`.
That license does not replace the licenses of source experimental data or
third-party software.

## Reinhardt RESI data

The Nup96 observations in `reproduce/data/witness_packing/observations.npz`
and the Nup96 display crop in `figures/figure_values.json` are reduced from:

Reinhardt SCM et al. *Ångström-resolution fluorescence microscopy* (2023),
[RESI data archive, doi:10.5281/zenodo.7795826](https://doi.org/10.5281/zenodo.7795826).
Source files are identified in `DATA.md` and the member index. The source
archive uses [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
The reduced files retain that license. Changes comprise selection of records,
coordinate and code extraction, group indexing, and construction of the witness
certificate and display crop. These files are not new experimental measurements.

## Other source records

Nup96 and CD20 sources are cited and linked in `DATA.md`. Their original
terms continue to apply to files obtained from those records.

## SPINNA and connected triple correlation

SPINNA is supplied by Picasso and is not redistributed as part of Crosslabel.
The reproduction scripts require `picassosr==0.11.3` and verify the upstream
module identity. The triple-correlation comparator is an independent scalar
implementation; its source article is cited in `DATA.md`.
