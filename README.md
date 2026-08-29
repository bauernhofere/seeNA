# ichorViz

<!-- badges: start -->
[![R-CMD-check](https://github.com/bauernhofere/ichorViz/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bauernhofere/ichorViz/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`ichorViz` is an independent R package for reading, comparing, and visualizing
existing [ichorCNA](https://github.com/GavinHaLab/ichorCNA) output. It adds a
file-based workflow for individual profiles, arbitrary genomic regions, paired
or longitudinal comparisons, and cohort-scale bin matrices and heatmaps.

It does **not** run ichorCNA or call copy-number alterations.

## Installation

```r
# install.packages("remotes")
remotes::install_github("bauernhofere/ichorViz")
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
plot_ichor_heatmap(mat, annotation_columns = c("condition", "timepoint"))
```

The included examples are generated format fixtures with no clinical or
biological interpretation.

## License and attribution

ichorViz is GPL-3-or-later and is not an official ichorCNA component. See
[NOTICE.md](NOTICE.md) for upstream provenance and attribution.
