# Data sources and reuse

The archive includes author-generated example catalogues, simulation
reference summaries, and derived analyses. Original third-party archives are available from their respective sources.
Reduced Nup96 observations used for the witness certificate and the figure
crop are included with attribution in THIRD_PARTY_NOTICES.md.

| Source | Public record | Role |
|---|---|---|
| Reinhardt et al., 2023 | [RESI archive](https://doi.org/10.5281/zenodo.7795826) | Nup96 measurements and CD20 provenance |
| Zaza et al., 2025 | [Nup96 archive](https://doi.org/10.6084/m9.figshare.28741373) | Second experimental Nup96 catalogue |
| Pachmayr et al., 2025 | [CD20 and antibody RESI archive](https://doi.org/10.5281/zenodo.15552356) | CD20 round-label provenance |
| Kowalewski et al., 2026 | [G5M archive](https://doi.org/10.5281/zenodo.18429674) | Public CD20 round-resolved source comparison |

Consult the source records for their terms. Crosslabel's MIT license applies
to its authored software and documentation. Original datasets and upstream
software retain the terms of their respective sources. SPINNA is part of
[Picasso](https://github.com/jungmannlab/picasso). The direct triple-correlation
comparator is an independent scalar adaptation of
[Yin et al.](https://doi.org/10.1038/s41467-018-08048-2).

## Nup96 source replay

Extract the `Figure2,EDF06,EDF07` directory from Reinhardt's `resi.rar`.
The replay uses `RESI/RESI rounds/RESI_R{1..4}_apicked.hdf5` and
`RESI/RESI result/RESI_R{1..4}_apicked_resi.hdf5` within that directory.
The [member index](reproduce/data/source_catalogues/cd20_members.json) lists the eight archive entries
under `npc`. For Zaza, extract `NUP96_resi.hdf5` from the Fig. 3 dataset.
The replay checks the Reinhardt array identities and the Zaza file hash
before scoring. See the [Nup96 command](reproduce/README.md#replay-the-experimental-catalogues).

## CD20 reconstruction comparison

From the G5M archive, retain the original four
`R{1..4}_picked_dbscan_molmap.hdf5` and
`R{1..4}_picked_ga_centers.hdf5` tables under
`data/Fig4/cell1_2023_Pos3/` and `data/Fig4/cell4_3D/`. The
[member index](reproduce/data/source_catalogues/g5m_members.json) records their exact offsets.
Check original member CRCs when extracting. Pass the extracted archive
root to the analysis. It verifies the scored G5M coordinates and the
gradient-ascent file hashes before evaluating the contrasts.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B \
  reproduce/field_only.py --cd20-sensitivity \
  --g5m-dir /path/to/extracted-g5m /path/to/cd20-sensitivity.json
```

## Included sensitivity results

- `reproduce/data/retention/`: all 2,400 fixed-footprint sensitivity fields,
  interval outputs, decisions, coverage, runtime measurements and study design.
- `reproduce/data/area/`: all 1,440 fixed-density area results, including
  480 reused baseline fields and 960 new fields. Catalogue hashes and seeds
  permit exact regeneration without storing duplicate raw catalogues.
- `reproduce/data/shared_interval/`: all 540 full-field triple-correlation
  interval results, sufficient arrays and nested calibration arrays.
- `reproduce/data/witness_packing/`: reduced Nup96 observations, group map
  and 110,732 record-disjoint witness triples. These observations derive from
  the Reinhardt RESI deposit; they retain CC BY 4.0 terms.

The manuscript figure inputs in `figures/figure_values.json` include a reduced
Nup96 display crop from the same deposit. Reproduction code and author-generated
simulation outputs are included in this release. Data from the other source
records remain available through the linked deposits above.
