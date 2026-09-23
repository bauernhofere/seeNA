.rebin_sample <- function(sample, bin_size, value, call_column) {
  d <- if (value == "segment_median") sample$segments else sample$bins
  column <- switch(value, call = call_column, segment_median = "median", value)
  if (!column %in% names(d)) .ichor_abort(paste("Required matrix column unavailable:", column))
  values <- if (value == "call") .call_score(d[[column]]) else d[[column]]
  d$value <- values
  d <- d[!is.na(d$value), c("chr", "start", "end", "value"), drop = FALSE]
  if (!nrow(d)) return(data.frame(chr = character(), bin_index = numeric(), value = numeric(), covered = numeric(), mixed = logical()))
  first <- floor((d$start - 1) / bin_size)
  last <- floor((d$end - 1) / bin_size)
  counts <- last - first + 1
  source_row <- rep.int(seq_len(nrow(d)), counts)
  bin_index <- if (all(counts == 1)) first else unlist(Map(seq.int, first, last), use.names = FALSE)
  overlap <- pmax(0, pmin(d$end[source_row], (bin_index + 1) * bin_size) -
                    pmax(d$start[source_row], bin_index * bin_size + 1) + 1)
  chr <- d$chr[source_row]
  v <- d$value[source_row]
  key <- paste(chr, bin_index, sep = ":")
  if (!anyDuplicated(key)) return(data.frame(chr = chr, bin_index = bin_index, value = v, covered = overlap, mixed = FALSE))
  keys <- unique(key)
  group <- match(key, keys)
  first_row <- match(keys, key)
  covered <- as.numeric(rowsum(overlap, group, reorder = FALSE))
  mixed <- rep(FALSE, length(keys))
  if (value == "call") {
    # Codes label categories, not a continuous scale. Aggregate base-pair support
    # by category; opposing gain/loss must never cancel into neutral.
    support <- vapply(c(-2, -1, 0, 1), function(code) {
      as.numeric(rowsum(overlap * (v == code), group, reorder = FALSE))
    }, numeric(length(keys)))
    support <- matrix(support, nrow = length(keys), ncol = 4L)
    max_support <- apply(support, 1L, max)
    tied <- rowSums(support == max_support) != 1L
    mixed <- rowSums(support > 0) > 1L
    aggregate_value <- c(-2, -1, 0, 1)[max.col(support, ties.method = "first")]
    aggregate_value[tied] <- NA_real_
  } else aggregate_value <- as.numeric(rowsum(overlap * v, group, reorder = FALSE)) / covered
  data.frame(chr = chr[first_row], bin_index = bin_index[first_row], value = aggregate_value,
             covered = covered, mixed = mixed)
}

#' Build a coverage-aware cohort matrix
#'
#' Rebins every sample in a cohort onto a fixed-width genome grid and returns
#' a sample-by-bin matrix of one chosen measurement together with matching
#' `coverage` and `mixed` layers.
#'
#' @details
#' Continuous values use base-pair-weighted means of the overlapping source
#' intervals. Calls use the category with the greatest base-pair support; ties give
#' `NA` and heterogeneous bins are flagged in `mixed`. Coverage is observed
#' overlap divided by target-bin width; cells below `min_coverage` or with no
#' observations are `NA`. Missing observations are never treated as neutral.
#'
#' `segment_median` reads the exported segment `median` directly, without changing
#' bin logR or applying a purity/ploidy transform. Every sample must have a segment
#' table. At boundaries, the result is an overlap-weighted mean of segment medians,
#' not a newly calculated median. NA medians and uncovered segment gaps contribute
#' no coverage. Coverage measures finite segment spans, not observed bin/read
#' coverage: a segment can bridge missing or filtered bins. No bin mask is applied.
#' These fitted summaries can smooth noise but do not validate low-TF fits or
#' make logR amplitudes directly comparable as tumor copy numbers across fluids.
#' See the installed methods contract for the formulas.
#' @param cohort An `ichor_cohort`.
#' @param bin_size Positive integer base-pair width.
#' @param value Explicit measurement: `"logR"`, `"corrected_copy_number"`,
#'   `"copy_number"`, `"call"`, or raw exported `"segment_median"` log2 ratio.
#' @param call_column For call matrices, `"corrected_call"` or `"event"`.
#' @param min_coverage Minimum observed fraction (default 1). Values with less
#'   support are `NA`; zero coverage always produces `NA`.
#' @param chromosomes Chromosomes to include, normalized and genomically sorted.
#' @param max_cells Allocation guard; maximum number of sample-by-bin cells.
#' @return A validated `ichor_matrix` with `values`, `coverage`, `mixed`,
#'   `bins`, `samples`, transformation settings and import provenance.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' cohort <- read_ichor_cohort(file.path(root, "example-manifest.csv"), "hg38")
#' m <- ichor_matrix(cohort, value = "call", chromosomes = "1")
#' m
#' m$values[, 1:4]
#' m$coverage[, 1:4]
#' logr <- ichor_matrix(cohort, value = "logR", bin_size = 5e5, chromosomes = "1")
#' logr$values[, 1:4]
#' segmented <- ichor_matrix(cohort, value = "segment_median", chromosomes = "1")
#' segmented$values[, 1:4]
#' @export
ichor_matrix <- function(cohort, bin_size = 1e6, value,
                         call_column = c("corrected_call", "event"), min_coverage = 1,
                         chromosomes = .chr_levels, max_cells = 5e7) {
  validate_ichor_cohort(cohort)
  value <- match.arg(value, c("logR", "corrected_copy_number", "copy_number", "call", "segment_median"))
  call_column <- match.arg(call_column)
  .scalar(bin_size, "bin_size", integer = TRUE, lower = 1)
  .scalar(min_coverage, "min_coverage", lower = 0, upper = 1)
  .scalar(max_cells, "max_cells", integer = TRUE, lower = 1)
  chromosomes <- unique(.normalize_chr(chromosomes))
  if (!length(chromosomes) || anyNA(chromosomes) || !all(chromosomes %in% .chr_levels)) .ichor_abort("Unsupported chromosomes.")
  source_layer <- if (value == "segment_median") "segments" else "bins"
  if (value == "segment_median") {
    missing <- vapply(cohort$samples, function(s) is.null(s$segments), logical(1))
    if (any(missing)) .ichor_abort(paste("Samples", paste(names(cohort$samples)[missing], collapse = ", "),
      "need segment tables for value = 'segment_median'; supply the selected-run segment files."), "seena_matrix_error")
  }
  sizes <- .chromosome_sizes(cohort$genome_build)
  sizes <- sizes[sizes$chr %in% chromosomes, ]
  n_bins <- sum(ceiling(sizes$length / bin_size))
  if (n_bins * length(cohort$samples) > max_cells) .ichor_abort("Matrix exceeds max_cells; increase bin_size or restrict chromosomes before allocation.")
  grid <- do.call(rbind, lapply(seq_len(nrow(sizes)), function(i) {
    idx <- seq_len(ceiling(sizes$length[i] / bin_size)) - 1
    data.frame(chr = sizes$chr[i], bin_index = idx, start = idx * bin_size + 1,
               end = pmin((idx + 1) * bin_size, sizes$length[i]))
  }))
  rownames(grid) <- NULL
  keys <- paste(grid$chr, grid$bin_index, sep = ":")
  grid$label <- sprintf("chr%s:%s-%s", grid$chr, grid$start, grid$end)
  dims <- list(names(cohort$samples), grid$label)
  mat <- matrix(NA_real_, length(cohort$samples), nrow(grid), dimnames = dims)
  coverage <- matrix(0, nrow(mat), ncol(mat), dimnames = dims)
  mixed <- matrix(FALSE, nrow(mat), ncol(mat), dimnames = dims)
  for (i in seq_along(cohort$samples)) {
    s <- cohort$samples[[i]]
    s[[source_layer]] <- s[[source_layer]][s[[source_layer]]$chr %in% chromosomes, , drop = FALSE]
    d <- .rebin_sample(s, bin_size, value, call_column)
    j <- match(paste(d$chr, d$bin_index, sep = ":"), keys)
    if (anyNA(j)) .ichor_abort("Source intervals are outside the target grid.")
    cov <- d$covered / (grid$end[j] - grid$start[j] + 1)
    d$value[cov + 1e-12 < min_coverage | cov == 0] <- NA_real_
    mat[i, j] <- d$value
    coverage[i, j] <- cov
    mixed[i, j] <- d$mixed
  }
  x <- structure(list(schema_version = 1L, values = mat, coverage = coverage, mixed = mixed,
                      bins = grid, samples = cohort$metadata, genome_build = cohort$genome_build,
                      bin_size = bin_size, value = value,
                      call_column = if (value == "call") call_column else NULL,
                      aggregation = if (value == "call") "base_pair_mode_ties_NA" else "base_pair_mean",
                      min_coverage = min_coverage, provenance = ichor_provenance(cohort),
                      coordinate_changes = lapply(cohort$samples, `[[`, "coordinate_changes")), class = "ichor_matrix")
  validate_ichor_matrix(x)
  x
}

#' Validate a cohort matrix
#'
#' Checks that values, coverage and mixed layers align with the bin grid and
#' sample metadata and obey the matrix measurement and coverage rules.
#' @param x An `ichor_matrix`.
#' @return `x`, invisibly.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' cohort <- read_ichor_cohort(file.path(root, "example-manifest.csv"), "hg38")
#' validate_ichor_matrix(ichor_matrix(cohort, value = "logR", chromosomes = "2"))
#' @export
validate_ichor_matrix <- function(x) {
  if (!inherits(x, "ichor_matrix") || !identical(x$schema_version, 1L) ||
      !is.matrix(x$values) || !is.numeric(x$values) || any(dim(x$values) == 0)) .ichor_abort("Invalid ichor_matrix.")
  if (!is.data.frame(x$bins) || !is.data.frame(x$samples) ||
      nrow(x$values) != nrow(x$samples) || ncol(x$values) != nrow(x$bins) ||
      !identical(rownames(x$values), x$samples$sample_id) || anyDuplicated(x$samples$sample_id) ||
      !identical(colnames(x$values), x$bins$label)) .ichor_abort("Matrix metadata or bin order disagrees.")
  .validate_intervals(x$bins, x$genome_build, c("chr", "start", "end"), "Matrix bins")
  .scalar(x$bin_size, "bin_size", integer = TRUE, lower = 1)
  if (length(x$value) != 1L || !x$value %in% c("call", "logR", "corrected_copy_number", "copy_number", "segment_median")) .ichor_abort("Invalid matrix measurement.")
  if (x$value == "call" && (length(x$call_column) != 1L || !x$call_column %in% c("corrected_call", "event"))) .ichor_abort("Invalid call column.")
  expected_aggregation <- if (x$value == "call") "base_pair_mode_ties_NA" else "base_pair_mean"
  if (!identical(x$aggregation, expected_aggregation)) .ichor_abort("Matrix aggregation disagrees with measurement.")
  grid <- x$bins
  sizes <- .chromosome_sizes(x$genome_build)
  limits <- sizes$length[match(grid$chr, sizes$chr)]
  if (!is.numeric(grid$bin_index) || anyNA(grid$bin_index) ||
      any(grid$bin_index != floor((grid$start - 1) / x$bin_size)) ||
      any(grid$start != grid$bin_index * x$bin_size + 1) ||
      any(grid$end != pmin((grid$bin_index + 1) * x$bin_size, limits)) ||
      !identical(grid$label, sprintf("chr%s:%s-%s", grid$chr, grid$start, grid$end))) .ichor_abort("Invalid target grid or labels.")
  .scalar(x$min_coverage, "min_coverage", lower = 0, upper = 1)
  for (nm in c("coverage", "mixed")) {
    if (!is.matrix(x[[nm]]) || !identical(dim(x[[nm]]), dim(x$values)) ||
        !identical(dimnames(x[[nm]]), dimnames(x$values))) .ichor_abort(paste("Invalid", nm, "matrix alignment."))
  }
  if (!is.numeric(x$coverage) || anyNA(x$coverage) || any(!is.finite(x$coverage) | x$coverage < 0 | x$coverage > 1 + 1e-12)) .ichor_abort("Invalid coverage fractions.")
  if (!is.logical(x$mixed) || anyNA(x$mixed)) .ichor_abort("Invalid mixed flags.")
  observed <- !is.na(x$values)
  if (any(!is.finite(x$values[observed])) ||
      any(observed & (x$coverage == 0 | x$coverage + 1e-12 < x$min_coverage))) .ichor_abort("Matrix values violate coverage or finiteness rules.")
  if (x$value == "call" && any(!x$values[observed] %in% c(-2, -1, 0, 1))) .ichor_abort("Fractional or unknown call codes.")
  if (x$value != "call" && any(x$mixed)) .ichor_abort("Mixed-call flags on a continuous matrix.")
  if (x$value %in% c("copy_number", "corrected_copy_number") && any(x$values[observed] < 0)) .ichor_abort("Negative matrix copy numbers.")
  invisible(x)
}

#' @export
print.ichor_matrix <- function(x, ...) {
  cat(sprintf("<ichor_matrix>\n  %d samples x %d bins\n  value: %s\n  bin size: %s bp\n  minimum coverage: %s\n",
              nrow(x$values), ncol(x$values), x$value, x$bin_size, x$min_coverage))
  invisible(x)
}
