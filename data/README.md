# Input data

The files in `data/raw/` are the inputs directly referenced by the standardized
MATLAB entry points. They are copied without numerical modification.

| File | Expected structure | Used by |
| --- | --- | --- |
| `phase_locking_inphase.xlsx` | Headerless numeric matrix: time, signal 1, signal 2 | `run_inphase_synchronization.m` |
| `phase_locking_antiphase.xlsx` | Headerless numeric matrix: time, signal 1, signal 2 | `run_antiphase_synchronization.m` |
| `ca1_in_vivo_recording.csv` | `Time`, `Intra_uV`, `Extra_uV` | `analyze_in_vivo_ca1_bursting.m` |
| `bnc_hardware_bursting.xlsx` | Headerless numeric matrix: time, hardware signal | `analyze_bnc_hardware_bursting.m` |
| `biological_segments.xlsx` | Time column followed by biological signal segments | `biomimetic_feature_overlap_analysis.m` |
| `hardware_segments.xlsx` | Time column followed by hardware signal segments | `biomimetic_feature_overlap_analysis.m` |

## Integrity

SHA-256 checksums are recorded in `SHA256SUMS.md`. Recompute them after any
intentional data replacement and review the resulting code outputs before
committing the change.

## Rights and redistribution

The data are included for the repository owner's research workflow. They are
not covered by the code license, and no permission for public redistribution is
granted by this repository. Keep the repository private until the owner has
confirmed the applicable source-dataset terms, participant protections, and
publication policy.

