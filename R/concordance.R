.pair_direction_class <- function(a, b) {
  out <- rep("unknown", length(a))
  ok <- !is.na(a) & !is.na(b)
  out[which(ok & a == 0 & b == 0)] <- "both_neutral"
  out[which(ok & a != 0 & b == 0)] <- "a_only"
  out[which(ok & a == 0 & b != 0)] <- "b_only"
  out[which(ok & a != 0 & b != 0 & sign(a) == sign(b))] <- "concordant"
  out[which(ok & a != 0 & b != 0 & sign(a) != sign(b))] <- "discordant"
  out
}

.pair_baseline_status <- function(s, call_column) {
  cn <- if (call_column == "corrected_call") "corrected_copy_number" else "copy_number"
  groups <- c("autosome", "X", "Y")
  if (!all(c(call_column, cn) %in% names(s$bins))) {
    return(stats::setNames(rep("unavailable_columns", 3), groups))
  }
  evidence <- ichor_neutral_cn(s, call_column)
  stats::setNames(evidence$status, evidence$chromosome_class)
}

.pair_joint_coverage <- function(a, b, call_column, grid, bin_size) {
  take <- function(s) {
    keep <- !is.na(.normalize_call(s$bins[[call_column]])) & s$bins$chr %in% grid$chr
    s$bins[keep, c("chr", "start", "end"), drop = FALSE]
  }
  a <- take(a); b <- take(b)
  n <- nrow(a) + nrow(b)
  chr <- character(n); start <- end <- numeric(n)
  ca <- match(a$chr, .chr_levels); cb <- match(b$chr, .chr_levels)
  i <- j <- 1L; k <- 0L
  # Sorted nonoverlapping inputs allow a linear intersection, not a Cartesian join.
  while (i <= nrow(a) && j <= nrow(b)) {
    if (ca[i] < cb[j]) { i <- i + 1L; next }
    if (cb[j] < ca[i]) { j <- j + 1L; next }
    left <- max(a$start[i], b$start[j]); right <- min(a$end[i], b$end[j])
    if (left <= right) {
      k <- k + 1L; chr[k] <- a$chr[i]; start[k] <- left; end[k] <- right
    }
    advance_a <- a$end[i] <= b$end[j]
    advance_b <- b$end[j] <= a$end[i]
    if (advance_a) i <- i + 1L
    if (advance_b) j <- j + 1L
  }
  covered <- numeric(nrow(grid))
  if (k) {
    rows <- seq_len(k)
    mask <- list(bins = data.frame(chr = chr[rows], start = start[rows], end = end[rows], logR = 0))
    d <- .rebin_sample(mask, bin_size, "logR", call_column)
    index <- match(paste(d$chr, d$bin_index), paste(grid$chr, grid$bin_index))
    covered[index] <- d$covered / (grid$end[index] - grid$start[index] + 1)
  }
  covered
}

.pair_representative <- function(d) {
  value <- rep(NA_real_, nrow(d)); source <- rep(NA_character_, nrow(d))
  put <- function(i, v, s) { value[i] <<- v; source[i] <<- s }
  put(which(d$concordance == "both_neutral"), 0, "none")
  for (s in c("a", "b")) {
    y <- d[[paste0("logR_", s)]]
    i <- which(d$concordance == paste0(s, "_only") & is.finite(y))
    put(i, y[i], s)
  }
  complete <- is.finite(d$logR_a) & is.finite(d$logR_b)
  i <- which(d$concordance == "concordant" & complete)
  put(i, (d$logR_a[i] + d$logR_b[i]) / 2, "mean")
  i <- which(d$concordance == "discordant" & complete)
  use_a <- abs(d$logR_a[i]) >= abs(d$logR_b[i])
  put(i, ifelse(use_a, d$logR_a[i], d$logR_b[i]), ifelse(use_a, "a", "b"))
  data.frame(representative_logR = value, representative_source = source)
}

#' Compare directional calls for an explicitly supplied pair
#'
#' Aligns two samples to the same reference grid using [ichor_matrix()]. Missing,
#' insufficiently covered, tied or mixed calls are unknown, not neutral. Calls
#' and logR have separate coverage: a known call can have an unavailable height.
#' Call support must also overlap over at least `min_coverage` of the target bin;
#' disjoint partially observed intervals are not compared as if co-observed.
#'
#' @details
#' Concordant means both altered in the same direction; loss and deep loss are
#' directionally concordant. Discordant means gain versus loss. This is not a
#' test of biological truth, magnitude agreement or fluid-specific detection.
#' Source calls are not reassigned based on logR sign or fitted TF.
#'
#' Sex-chromosome baseline evidence is always reported. Default `"flag"` retains
#' supported source calls even if no NEUT bins exist: a genuinely altered whole
#' chromosome can lack NEUT bins. `"require_neutral"` instead makes X/Y comparisons
#' unknown when either sample lacks a unique observed NEUT-bin CN. Neither policy
#' verifies a biological baseline; no autosomal fallback or assumed CN is used.
#'
#' `representative_logR` reproduces the manuscript display rule: the altered
#' sample for one-sided calls, mean for concordant calls, and larger absolute logR
#' for discordant calls (ties choose a). Both heights must exist for a two-sided
#' representative. This quantity is not a difference or effect size. Swapping
#' samples can reverse its sign in equal-absolute discordant ties.
#' @param a,b Explicit `ichor_sample` objects, with distinct IDs and the same build.
#' @param bin_size Target bin width in bp; default 1 Mb.
#' @param call_column `"corrected_call"` or `"event"`; no fallback.
#' @param min_coverage Required observed fraction for each call/logR layer and
#'   joint call support. Partial logR means can still summarize different observed
#'   portions; heights are descriptive summaries, not a paired difference.
#' @param chromosomes Reference chromosomes to include, default 1-22/X.
#' @param ploidy_adjust Apply the upstream shift to aligned logR? Default FALSE.
#' @param sex_chromosomes `"flag"` or conservative `"require_neutral"` policy.
#' @param max_cells Allocation guard passed to each two-row matrix builder.
#' @return Data frame in genomic order, including coordinates, aggregated display
#'   calls, raw/selected-scale logR, separate coverage, mixed flags, baseline status,
#'   concordance/reason, and representative height/source. All target bins remain,
#'   including both-neutral and unknown bins. Attributes `settings`, `provenance`
#'   and `coordinate_changes` retain the import and transformation context.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' a <- read_ichor_sample(file.path(root, "example-a.cna.seg"), genome_build = "hg38")
#' b <- read_ichor_sample(file.path(root, "example-b.cna.seg"), genome_build = "hg38")
#' d <- ichor_pair_concordance(a, b, chromosomes = "1")
#' d[1:4, c("chr", "start", "call_a", "call_b", "concordance")]
#' @export
ichor_pair_concordance <- function(a, b, bin_size = 1e6,
                                    call_column = c("corrected_call", "event"),
                                    min_coverage = 1,
                                    chromosomes = c(as.character(1:22), "X"),
                                    ploidy_adjust = FALSE,
                                    sex_chromosomes = c("flag", "require_neutral"),
                                    max_cells = 5e7) {
  samples <- .as_sample_list(list(a, b))
  if (!identical(a$genome_build, b$genome_build)) .ichor_abort("Pair must use the same genome build.")
  call_column <- match.arg(call_column)
  sex_chromosomes <- match.arg(sex_chromosomes)
  .flag(ploidy_adjust, "ploidy_adjust")
  shifts <- if (ploidy_adjust) vapply(samples, .logr_shift, numeric(1)) else c(0, 0)
  cohort <- structure(list(schema_version = 1L, samples = samples,
    metadata = data.frame(sample_id = names(samples)), genome_build = a$genome_build), class = "ichor_cohort")
  calls <- ichor_matrix(cohort, bin_size, value = "call", call_column = call_column,
    min_coverage = min_coverage, chromosomes = chromosomes, max_cells = max_cells)
  logr <- ichor_matrix(cohort, bin_size, value = "logR", min_coverage = min_coverage,
    chromosomes = chromosomes, max_cells = max_cells)
  d <- calls$bins
  state_names <- c("-2" = "Deep loss", "-1" = "Loss", "0" = "Neutral", "1" = "Gain")
  chr_class <- ifelse(d$chr %in% c("X", "Y"), d$chr, "autosome")
  for (i in 1:2) {
    s <- c("a", "b")[i]
    d[[paste0("call_", s)]] <- unname(state_names[as.character(calls$values[i, ])])
    d[[paste0("coverage_", s)]] <- unname(calls$coverage[i, ])
    d[[paste0("mixed_", s)]] <- unname(calls$mixed[i, ])
    d[[paste0("raw_logR_", s)]] <- unname(logr$values[i, ])
    d[[paste0("logR_", s)]] <- unname(logr$values[i, ]) + shifts[i]
    d[[paste0("logR_coverage_", s)]] <- unname(logr$coverage[i, ])
    d[[paste0("baseline_status_", s)]] <- unname(.pair_baseline_status(samples[[i]], call_column)[chr_class])
  }
  d$joint_coverage <- .pair_joint_coverage(a, b, call_column, d, bin_size)
  d$baseline_flag <- d$chr %in% c("X", "Y") &
    (d$baseline_status_a != "supported" | d$baseline_status_b != "supported")
  d$concordance <- .pair_direction_class(calls$values[1, ], calls$values[2, ])
  d$reason <- rep("", nrow(d))
  add_reason <- function(mask, text) {
    i <- which(mask)
    d$reason[i] <<- ifelse(nzchar(d$reason[i]), paste(d$reason[i], text, sep = ";"), text)
  }
  add_reason(is.na(d$call_a), "call_unavailable_a")
  add_reason(is.na(d$call_b), "call_unavailable_b")
  add_reason(d$mixed_a, "mixed_calls_a")
  add_reason(d$mixed_b, "mixed_calls_b")
  add_reason(d$joint_coverage == 0 | d$joint_coverage + 1e-12 < min_coverage, "insufficient_joint_call_support")
  if (sex_chromosomes == "require_neutral") add_reason(d$baseline_flag, "sex_reference_unresolved")
  d$concordance[nzchar(d$reason)] <- "unknown"
  d <- cbind(d, .pair_representative(d))
  attr(d, "settings") <- list(schema_version = 1L, sample_ids = names(samples),
    genome_build = a$genome_build, bin_size = bin_size, call_column = call_column,
    min_coverage = min_coverage, chromosomes = unique(d$chr),
    ploidy_adjust = ploidy_adjust, shifts = stats::setNames(shifts, names(samples)),
    mixed_calls = "unknown", sex_chromosomes = sex_chromosomes,
    representative = "one_sided_or_mean_or_largest_absolute_ties_a")
  attr(d, "provenance") <- calls$provenance
  attr(d, "coordinate_changes") <- calls$coordinate_changes
  d
}
