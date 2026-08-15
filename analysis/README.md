# Analysis code

## Standardized filename mapping

| Source filename | Repository filename |
| --- | --- |
| `suoxinag.m` | `run_antiphase_synchronization.m` |
| `suoxinag2.m` | `run_inphase_synchronization.m` |
| `Biomimetic_Indistinguishability_Analysis.m` | `biomimetic_feature_overlap_analysis.m` |
| `InVivo_CA1_Bursting_Analyzer.m` | `analyze_in_vivo_ca1_bursting.m` |
| `Neuromorphic_HW_Bursting_Analyzer.m` | `analyze_bnc_hardware_bursting.m` |

The standardized names describe the implemented operation without making a
scientific interpretation part of the filename.

## Entry points

Run scripts from the repository root with MATLAB's `run` function. Each entry
point constructs repository-relative input and output paths, then calls a
reusable analysis function where appropriate.

## Shared functions

- `compute_phase_frequency_map.m` implements filtering, Hilbert-domain phase
  extraction, phase-binned cross-power, coherence estimation, visualization,
  and structured output.
- `run_burst_analysis.m` implements robust normalization, polarity-aware spike
  detection, burst grouping, ISI metrics, return maps, and tabular exports.

## Output contract

Scripts write only beneath `results/`. They do not overwrite files in
`data/raw/` and do not depend on the original local source-code directory.

