# Agent guide

## Scope

ichorViz reads and visualizes existing ichorCNA outputs. It does not call CNAs,
select an ichorCNA solution, or import clinical spreadsheets.

## Architecture

`files -> parsers -> ichor_sample -> ichor_cohort -> aligned matrix -> plots`

- `R/read.R`: file contracts and parsers
- `R/sample.R`: sample object and validation
- `R/genome.R`: builds, coordinates, and regions
- `R/plot.R`: individual and comparison plots
- `R/cohort.R`: manifest import and harmonized matrices
- `R/heatmap.R`: optional ComplexHeatmap adapter

## Invariants

- Inputs are read-only.
- Internal tumor fractions are fractions in `[0, 1]`.
- Genome build is explicit and is never guessed.
- Chromosomes are normalized to `1`-`22`, `X`, and `Y`.
- Parsing, scientific transformations, and rendering stay separate.
- Plot functions return objects; they do not write files.
- No patient-specific corrections, absolute paths, or clinical data.
- Tests do not require network access.
- Never commit real `.cna.seg`, `.seg`, or `.params.txt` files.

## Verification

Run `make check`. It generates documentation, runs tests, and performs
`R CMD check`. Changes are complete only when this command succeeds without
warnings.
