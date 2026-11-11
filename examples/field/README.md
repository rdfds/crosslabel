# Example catalogues

These synthetic catalogues are the first two cases of the prespecified
held-out matched-moment study. Their expected JSON reports are checked from
an installed wheel. Use one numerical thread when comparing saved values.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 crosslabel-field \
  --csv examples/field/matched_adequate.csv --radius-nm 30 --seed 984000
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 crosslabel-field \
  --csv examples/field/matched_alternative.csv --radius-nm 30 --seed 984500
```

Compare the estimated values and interval limits with the corresponding
JSON files. The saved output was generated on Python 3.13 with NumPy and
SciPy. Minor rounding differences can occur across numerical environments.
No public experimental coordinate data are included in these examples.
