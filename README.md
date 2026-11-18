# Crosslabel

Crosslabel checks whether a model built from single molecules and pairs can
explain how three imaging codes occur together in a localization microscopy
field. It measures the three-code signal and estimates its uncertainty.

## Install

Use Python 3.9 or later. From the repository root, install the bundled wheel
with HDF5 support.

```sh
python -m pip install './dist/crosslabel-0.4.1-py3-none-any.whl[hdf5]'
```

## Quickstart

Run the field check on a synthetic catalogue.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 crosslabel-field \
  --csv examples/field/matched_alternative.csv --radius-nm 30 --seed 984500
```

The same check is available in Python.

```python
from crosslabel import field_check
from crosslabel.io import read_table

xy, codes, names = read_table(
    "examples/field/matched_alternative.csv", x="x_nm", y="y_nm", label="code"
)
result = field_check(xy, codes, code_count=len(names), radius_nm=30, seed=984500)
print(result["estimate"], result["interval_95_model"])
```

## Inputs and outputs

Use a CSV with `x_nm,y_nm,code`, or one Picasso HDF5 file per imaging code
with a pixel size in nanometres. The validated field configurations contain **4 or 6
permanent codes**. Choose a neighbourhood radius of **10 or 30 nm**.

The report gives the contrast, a model-based interval, and checks for
sampling support and spatial density. See the [field guide](reproduce/FIELD_TOOL.md)
for the options and [example outputs](examples/field/).

## Interpretation

The interval describes uncertainty under the stated component and observation
models. An experimental decision requires evidence that the interval or a
simulation reference bank is calibrated for the acquisition.

## Reproduce and cite

The [reproduction guide](reproduce/README.md) links to the results and gives
commands for the checks, replays, and figures. [Data sources](DATA.md) lists the
original deposits. Simulation results and reference tables are included in
[`reproduce/data/`](reproduce/data/).

See [CITATION.cff](CITATION.cff) for the software citation. Crosslabel is
released under the [MIT license](LICENSE).

Version 0.4.1, including all reported sensitivity studies, is archived at
[doi:10.5281/zenodo.24354693](https://doi.org/10.5281/zenodo.24354693).
