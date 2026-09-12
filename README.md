# ichorViz

<!-- badges: start -->
[![R-CMD-check](https://github.com/bauernhofere/ichorViz/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bauernhofere/ichorViz/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`ichorViz` is an independent R package for reading, comparing, and visualizing
existing [ichorCNA](https://github.com/GavinHaLab/ichorCNA) output. It adds a
file-based workflow for individual profiles, arbitrary genomic regions, paired
or longitudinal comparisons, and cohort-scale bin matrices and heatmaps.

It does **not** run ichorCNA or call copy-number alterations.

**Pre-release:** the repository is private and this development version is not
an archived/citable release. See [scientific contracts](docs/data-contracts.md)
and [review/release gates](docs/review-triage.md). The checks below do not replace
real-input reconciliation and author review of manuscript figures.

## Installation

```r
# install.packages("remotes")
# Requires access to the private repository; pin an approved commit for analysis.
remotes::install_github("bauernhofere/ichorViz")

# Optional cohort heatmaps:
# install.packages("BiocManager")
# BiocManager::install(c("ComplexHeatmap", "circlize"))
# install.packages("ragg")
```

## One sample

```r
library(ichorViz)

root <- system.file("extdata", package = "ichorViz")
sample_a <- read_ichor_sample(
  file.path(root, "example-a.cna.seg"),
  file.path(root, "example-a.seg"),
  file.path(root, "example-a.params.txt"),
  genome_build = "hg38"
)

plot_ichor_profile(sample_a)
plot_ichor_region(sample_a, "chr1:1-4000000")
```

## Compare related samples

```r
sample_b <- read_ichor_sample(
  file.path(root, "example-b.cna.seg"),
  file.path(root, "example-b.seg"),
  file.path(root, "example-b.params.txt"),
  genome_build = "hg38"
)

plot_ichor_compare(
  list(sample_a, sample_b),
  colors = c("example-a" = "#8E0D0D", "example-b" = "#EDB332")
)

plot_ichor_region(
  list(sample_a, sample_b),
  region = "chrX:1-4000000",
  colors = c("example-a" = "#8E0D0D", "example-b" = "#EDB332")
)
```

## Cohort heatmap

```r
cohort <- read_ichor_cohort(
  file.path(root, "example-manifest.csv"),
  genome_build = "hg38"
)

mat <- ichor_matrix(cohort, bin_size = 1e6, value = "corrected_copy_number")
h <- plot_ichor_heatmap(mat, annotation_columns = c("condition", "timepoint"))
ComplexHeatmap::draw(h)
```

The included examples are generated format fixtures with no clinical or
biological interpretation.

## Scientific defaults

- `value` is required for a matrix; `call_column` never falls back silently.
- Missing calls are NA; unsupported tokens are errors, not neutral calls.
- Continuous values use overlap-weighted means. Categorical calls use base-pair
  mode (ties -> NA); inspect `mat$mixed` for heterogeneous bins.
- `mat$coverage` reports observed support. Default `min_coverage = 1` requires
  full target coverage, including shortened chromosome-end bins.
- Genome bounds are strict. If you have verified terminal bin padding in the
  source output, use `bounds = "trim"` at import and inspect `coordinate_changes`.
- Heatmaps retain manifest order; clustering is explicit and uses common observed
  bins without imputation. CN 2 is a visual reference, not an inferred baseline.
- Source component IDs must match before aliases are applied. You must choose
  the same ichorCNA run for all components; matching IDs cannot verify the run.
- Provenance fingerprints are captured at import; paths require
  `retain_paths = TRUE`. Neither aliases nor path redaction anonymize genomic data.

## Development and manuscript use

```bash
make bootstrap  # installs development dependencies; network required
make test
make check
```

Use real clinical inputs outside this repository. Record `sessionInfo()`, the
package commit, `ichor_provenance(cohort)`, the selected run, genome build,
coordinate changes, matrix settings and plotting options alongside manuscript
outputs in an access-controlled location. No public data upload is required.
For an accessible software citation, an author-approved public release/archive
is still needed; cite the original ichorCNA paper separately:
[Adalsteinsson, Ha, Freeman et al. (2017)](https://doi.org/10.1038/s41467-017-00965-y).

## License and attribution

ichorViz is GPL-3-or-later and is not an official ichorCNA component. See
[NOTICE.md](NOTICE.md) for upstream provenance and attribution. A copy ships
with the installed package at `system.file("NOTICE.md", package = "ichorViz")`.
