# Reproduce the results

Install the wheel with its HDF5 extra and run these commands from the
repository root. The forward replays were checked on Python 3.13 with
NumPy 2.5.3, SciPy 1.18.1, and h5py 3.16.0. Use one thread for numerical
comparisons. To install the original numerical dependencies in a fresh
Python 3.13 environment, run the following from the repository root.

```sh
python -m pip install 'numpy==2.5.3' 'scipy==1.18.1' 'h5py==3.16.0' \
  './dist/crosslabel-0.4.0-py3-none-any.whl'
```

## Verify the tables

The matched-moment replay regenerates the pair training bank, the reference
bank, and both methods' fits and decisions on the held-out fields. The other
commands recompute decisions and summaries from saved scores. Binomial
interval endpoints use a relative tolerance of `1e-9`. Counts, ranks, and
input hashes require exact agreement.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -B reproduce/replay_marked_pairs.py \
  reproduce/matched_pair_calibration.npz reproduce/data/marked_pair/MATCHED.jsonl
python -B reproduce/check_robustness.py reproduce/data/robustness --check
python -B reproduce/match_false_alerts.py reproduce/data/robustness | \
  cmp - reproduce/data/robustness/OPERATING_POINTS.csv
python -B reproduce/check_nup96.py reproduce/data/nup96 --check
```

## Replay the simulations

Use new output paths. The robustness replay regenerates the pair training
fields, the reference and evaluation fields, and the Crosslabel,
triple-correlation, and marked-pair scores. It checks the regenerated scores
and decisions against the corresponding saved results. The field-only
replay also computes Crosslabel intervals for the uniform intensity controls.

```sh
python -B reproduce/collect_robustness.py /path/to/new-robustness
python -B reproduce/field_only.py /path/to/field-only.json \
  --robustness-dir /path/to/new-robustness
python -B reproduce/tile_sensitivity.py /path/to/tile-sensitivity.json \
  --field-results /path/to/field-only.json
```

The interval replay regenerates 400 native fields in each of six observation
conditions and 20 paired merge controls. It computes both code masks and
all nine interval variants and verifies the coverage and paired comparison
tables. Use multiple workers if available, with one numerical thread per
worker. The full run can take substantial time. Add `--smoke` for an
eight-view installation check.

```sh
python -B reproduce/replay_intervals.py /path/to/new-intervals --workers 4
```

## Replay the experimental catalogues

The [Nup96 input map](../DATA.md#nup96-source-replay) identifies the original
HDF5 files. The command recomputes scores and intervals for both native
catalogues, 99 reference fields, and 200 paired template fields, then compares
them to the saved values and ranks.

```sh
python -B reproduce/replay_nup96.py \
  --reinhardt '/path/to/Figure2,EDF06,EDF07' --zaza /path/to/NUP96_resi.hdf5
```

For the CD20 point contrasts and subsampling, use the
[CD20 input map](../DATA.md#cd20-reconstruction-comparison).

## SPINNA comparison

The [SPINNA results](data/spinna/REPORT.md) include the reported decisions and
fitted-null comparisons. The [method](spinna_comparison.py) uses
`picassosr==0.11.3` and checks the hash of the upstream SPINNA module.
Install that version of Picasso into the replay environment. The commands
regenerate the shared inputs, rerun both methods, and check the saved
paired decisions and SPINNA ranks. Each example runs one field. Omit
`--seed` for all 144 or 160 paired fields, respectively. The full refits
are substantially slower than the table checks.

```sh
python -m pip install 'picassosr==0.11.3'
python -B reproduce/replay_spinna.py geometry /path/to/new-geometry --seed 801001
python -B reproduce/replay_spinna.py library /path/to/new-library --seed 701001
```

To verify both masks and all interval variants on the 1,200 geometry
fields or 160 library fields, use `--intervals-only` and omit `--seed`.
This mode does not require Picasso.

## Generate the figures

Install the wheel with the `figures` extra and supply a new output directory.

```sh
python -B reproduce/plot_workflow.py /path/to/new-figures
python -B reproduce/plot_benchmarks.py /path/to/new-figures/viability_rates.pdf
```

## Results

- [Marked-pair benchmark and Nup96 controls](data/marked_pair/REPORT.md)
- [Interval qualification](data/interval_qualification/REPORT.md)
- [Robustness and triple-correlation comparison](data/robustness/REPORT.md)
- [Matched observed false-alert analysis](data/robustness/MATCHED_FALSE_ALERTS.md)
- [Nup96 comparison](data/nup96/REPORT.md)
- [Field-only intervals, intensity control, and CD20 sensitivity](data/field_only/REPORT.md)
- [Experimental acquisition audit](data/acquisition_audit/REPORT.md)
