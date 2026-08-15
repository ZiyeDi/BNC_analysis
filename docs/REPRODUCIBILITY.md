# Reproducibility checklist

## Environment

Record the following before generating release results:

- MATLAB release and operating system
- Signal Processing Toolbox version
- Statistics and Machine Learning Toolbox version
- Repository commit identifier
- Input-data checksums

## Execution

1. Confirm that all files listed in `data/SHA256SUMS.md` match their recorded
   checksums.
2. Start MATLAB in the repository root.
3. Run one analysis entry point at a time.
4. Preserve the complete console log for the run.
5. Check the generated PNG and FIG files for missing panels or labels.
6. Confirm that each CSV column has the documented unit and that MAT files load
   without warnings.

## Interpretation boundary

The software reports numerical outputs for configured datasets and parameters.
Repository documentation does not make manuscript-level claims or substitute
for statistical review. Any public interpretation should be written and
approved separately from this private code repository.

## Release preparation

Before changing repository visibility or creating a release:

- confirm data redistribution rights;
- remove any file that is not cleared for release;
- update `CITATION.cff` with the final version and publication metadata;
- tag the validated commit; and
- archive the exact MATLAB environment and generated outputs separately.

