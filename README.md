# BNC_analysis

Private MATLAB research-code repository for paired-signal synchronization,
bursting analysis, and selected energy-envelope feature comparisons.

This repository intentionally excludes the manuscript, manuscript excerpts,
figures, unpublished narrative, and reported result values. Its documentation
is limited to code interfaces, input schemas, execution steps, and output
artifacts.

## Analysis modules

| Script | Purpose |
| --- | --- |
| `analysis/run_inphase_synchronization.m` | Phase-frequency mapping for the in-phase paired-signal dataset. |
| `analysis/run_antiphase_synchronization.m` | Phase-frequency mapping for the anti-phase paired-signal dataset. |
| `analysis/analyze_in_vivo_ca1_bursting.m` | Spike and burst analysis for the CA1 extracellular recording. |
| `analysis/analyze_bnc_hardware_bursting.m` | Spike and burst analysis for the BNC hardware signal. |
| `analysis/biomimetic_feature_overlap_analysis.m` | Selected energy-envelope features and cross-validated linear-SVM evaluation. |

Shared functions are kept in `analysis/compute_phase_frequency_map.m` and
`analysis/run_burst_analysis.m` to reduce duplicated logic.

## Repository layout

```text
BNC_analysis/
|-- analysis/          MATLAB entry points and reusable functions
|-- data/raw/          Input files referenced by the entry-point scripts
|-- docs/              Reproducibility and maintenance notes
|-- results/           Generated outputs (not versioned by default)
|-- CITATION.cff       Repository citation metadata
|-- CONTRIBUTING.md    Contribution workflow
|-- LICENSE            Code and data rights notice
`-- README.md
```

## Requirements

- MATLAB with `readmatrix`, `writetable`, and `savefig`
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

The code uses only repository-relative paths. No source file refers to the
original local drive locations.

## Quick start

1. Clone the repository and open MATLAB in the repository root.
2. Run one entry-point script, for example:

   ```matlab
   run(fullfile('analysis', 'analyze_in_vivo_ca1_bursting.m'))
   ```

3. Review generated figures, tables, reports, and MAT files under `results/`.

To run another workflow, replace the entry-point filename with one from the
analysis-module table above.

## Reproducibility notes

- Raw input files are treated as immutable.
- Each script resolves inputs from `data/raw/` and writes outputs to a dedicated
  folder under `results/`.
- Generated outputs are ignored by Git to prevent accidental commits of local
  reruns.
- Input file checksums are recorded in `data/SHA256SUMS.md`.

See `docs/REPRODUCIBILITY.md` for the full execution and verification checklist.

## Data and confidentiality

The repository is intended to remain private unless the owner completes a
separate publication and data-rights review. The manuscript is not part of this
repository. Data files are not relicensed by the code license; see
`data/README.md` and `LICENSE`.

## Citation

Use the metadata in `CITATION.cff` when citing this code repository. Update the
version, release date, and publication identifiers before a public release.

