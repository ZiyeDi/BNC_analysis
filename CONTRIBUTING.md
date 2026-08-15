# Contributing

This is a private research repository. Coordinate changes with the repository
owner before opening a pull request.

## Change requirements

- Keep all code comments, messages, labels, and filenames in English.
- Use repository-relative paths; do not commit machine-specific drive paths.
- Do not modify files under `data/raw/` without updating `data/SHA256SUMS.md`.
- Do not add manuscript text, manuscript figures, or unpublished narrative.
- Run the affected MATLAB entry point and inspect all generated figures.
- Summarize parameter or algorithm changes in the pull-request description.

## Commit scope

Prefer focused commits that separate code changes, data replacements, and
documentation updates. Never commit generated `results/` files unless the
repository owner explicitly approves a release artifact.

