# Crosslabel field check

Crosslabel estimates a connected three-code contrast from a full catalogue of
2D molecular localizations or molecule positions. It returns the model interval
used in the manuscript and a separate spatial density diagnostic. It takes no
assembly or pore assignments. Its inputs are positions and permanent codes.

Install from this source tree in a Python 3.9 or later environment.

```sh
python -m pip install .
crosslabel-field --csv examples/field/matched_alternative.csv --radius-nm 30 --seed 984500
```

The CSV header is `x_nm,y_nm,code`. Use one code per imaging round and input
all rows in the declared catalogue. The command uses either four or six codes,
with 10 or 30 nm neighbourhoods. It records its numeric seed and window.

To use one Picasso HDF5 file per code, install the HDF5 extra and specify the
camera pixel size independently of the localization data.

```sh
python -m pip install '.[hdf5]'
crosslabel-field --picasso round1.hdf5 round2.hdf5 round3.hdf5 round4.hdf5 \
  --pixel-nm 130 --radius-nm 30 --seed 42 --output field.json
```

Python users can call `crosslabel.field.field_check` on a NumPy array of
`[x_nm,y_nm]` rows and an integer code array. Supply `code_count`,
`radius_nm`, and `seed` as keywords. The output has the estimate, interval,
numerical support, and density score. A positive lower limit indicates
incompatibility under the stipulated independent-component model and the
stated observation assumptions. Separate evidence for the acquisition and
calibration is required before using this value to reject an experimental
library. A nonpositive lower limit leaves the library unresolved.

For the density diagnostic, use `crosslabel.field_screen.density_score(xy)`
to collect one score per independent reference field. Reference fields should
match the support, density, sampling and observation process expected under
the proposed model. Pass their scores as a JSON array using
`--reference-scores`. A reference rank at or below .05 flags excess coarse
variation for review. Passing this screen provides no coverage guarantee.
Check the common footprint, permanent code assignments and preprocessing
against acquisition records and relevant calibration controls.

The supplied examples are synthetic controls. The two fixed seeds belong to
the matched-moment follow-up in the manuscript. Expected JSONs and comparison
commands are in `examples/field/`. The CLI and Python API use the same field estimator.

The code is MIT licensed. Cite Reinhardt et al. for the Nup96 source data
([Nature article](https://doi.org/10.1038/s41586-023-05925-9) and
[dataset](https://doi.org/10.5281/zenodo.7795826)).
