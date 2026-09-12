#' Read the fitted tumor fraction
#' @param x An ichor_sample.
#' @return Numeric scalar between zero and one, or NA_real_ if unavailable.
#' @export
ichor_tf <- function(x) {
  validate_ichor_sample(x)
  if (is.null(x$params)) NA_real_ else x$params$tumor_fraction
}

#' Read fitted tumor ploidy
#' @param x An ichor_sample.
#' @return Positive numeric scalar, or NA_real_ if parameters are unavailable.
#' @export
ichor_ploidy <- function(x) {
  validate_ichor_sample(x)
  if (is.null(x$params)) NA_real_ else x$params$ploidy
}

.logr_shift <- function(x) {
  p <- x$params
  if (is.null(p) || !is.finite(p$tumor_fraction) || !is.finite(p$ploidy)) {
    .ichor_abort("Ploidy-adjusted logR requires finite tumor fraction and ploidy from the selected run.", "ichorviz_parameter_error")
  }
  log2((p$tumor_fraction * p$ploidy + (1 - p$tumor_fraction) * 2) / 2)
}

#' Apply the upstream ichorCNA plotting shift
#'
#' Adds `log2((TF * ploidy + (1 - TF) * 2) / 2)` to source logR or segment
#' medians, matching the formula in v0.4 plotGWSolution/plotCNlogRByChr. Exported
#' parameters are rounded, so results need not be identical to a plot made from
#' full-precision in-memory fitted parameters. This is a display transform, not
#' purity correction, integer copy-number estimation or call reassignment. It
#' does not guarantee that every gain/loss lies above/below zero, particularly
#' on sex chromosomes or in noisy bins. Raw input values are never modified.
#' @param x An ichor_sample with finite fitted TF and ploidy.
#' @param component Bins (logR) or segments (median).
#' @return Numeric vector in source-table order, retaining NA. For an absent
#'   segment table, returns numeric(0); fitted parameters are still required.
#' @export
ichor_adjusted_logr <- function(x, component = c("bins", "segments")) {
  validate_ichor_sample(x)
  component <- match.arg(component)
  shift <- .logr_shift(x)
  values <- if (component == "bins") x$bins$logR else x$segments$median
  as.numeric(values) + shift
}

#' Summarize copy numbers observed in neutral-called bins
#'
#' Reports empirical evidence, not a guaranteed biological or run-wide neutral
#' baseline. Exactly one distinct observed CN within a chromosome class yields
#' a value; absent evidence or multiple distinct values yields NA with a status.
#' Never defaults to two copies, uses gender to guess a baseline, or borrows an
#' autosomal value for X/Y. Call and CN columns are paired from the same layer.
#' @param x An ichor_sample.
#' @param call_column corrected_call (paired with corrected_copy_number) or event
#'   (paired with copy_number).
#' @return Three-row data frame for autosome, X and Y: chromosome_class,
#'   neutral_cn, status, n_neutral_bins, n_observed, n_distinct, and a list-column
#'   candidates. Status is supported, no_neutral_bins, missing_cn, or ambiguous.
#'   Check status before using the result in a downstream transform.
#' @export
ichor_neutral_cn <- function(x, call_column = c("corrected_call", "event")) {
  validate_ichor_sample(x)
  call_column <- match.arg(call_column)
  cn_column <- if (call_column == "corrected_call") "corrected_copy_number" else "copy_number"
  d <- x$bins
  if (!all(c(call_column, cn_column) %in% names(d))) .ichor_abort("Neutral-CN evidence requires the paired call and CN columns.")
  chr_class <- ifelse(d$chr %in% as.character(1:22), "autosome", d$chr)
  calls <- .normalize_call(d[[call_column]])
  out <- lapply(c("autosome", "X", "Y"), function(group) {
    values <- d[[cn_column]][which(chr_class == group & calls == "NEUT")]
    candidates <- sort(unique(values[!is.na(values)]))
    status <- if (!length(values)) "no_neutral_bins" else if (!length(candidates)) "missing_cn" else
      if (length(candidates) > 1) "ambiguous" else "supported"
    result <- data.frame(chromosome_class = group,
      neutral_cn = if (status == "supported") candidates else NA_real_, status = status,
      n_neutral_bins = length(values), n_observed = sum(!is.na(values)), n_distinct = length(candidates))
    result$candidates <- list(candidates)
    result
  })
  do.call(rbind, out)
}
