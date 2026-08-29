# Sample object ------------------------------------------------------------

#' Construct a validated ichorCNA sample
#'
#' @param cna_seg Path to a `.cna.seg` file.
#' @param seg Optional path to a `.seg` file.
#' @param params Optional path to a `.params.txt` file.
#' @param genome_build Genome build, currently `"hg19"` or `"hg38"`.
#' @param sample_id Optional sample identifier.
#' @return An object of class `ichor_sample`.
#' @export
read_ichor_sample <- function(cna_seg, seg = NULL, params = NULL,
                              genome_build, sample_id = NULL) {
  genome_build <- match.arg(genome_build, c("hg19", "hg38"))
  bins <- read_ichor_cna(cna_seg, sample_id = sample_id)
  segments <- if (is.null(seg) || is.na(seg) || !nzchar(seg)) NULL else read_ichor_segments(seg)
  parameters <- if (is.null(params) || is.na(params) || !nzchar(params)) NULL else read_ichor_params(params)
  id <- sample_id %||% attr(bins, "sample_id") %||%
    if (!is.null(parameters)) parameters$sample_id[1] else NULL

  x <- structure(
    list(
      sample_id = as.character(id),
      genome_build = genome_build,
      bins = bins,
      segments = segments,
      params = parameters,
      provenance = list(
        cna_seg = attr(bins, "source"),
        seg = if (!is.null(segments)) attr(segments, "source") else NULL,
        params = if (!is.null(parameters)) normalizePath(params, mustWork = TRUE) else NULL
      )
    ),
    class = "ichor_sample"
  )
  validate_ichor_sample(x)
  x
}

#' Validate an ichorCNA sample
#'
#' @param x An `ichor_sample` object.
#' @return `x`, invisibly. Invalid samples raise an error.
#' @export
validate_ichor_sample <- function(x) {
  if (!inherits(x, "ichor_sample")) .ichor_abort("x must be an ichor_sample.")
  if (length(x$sample_id) != 1L || is.na(x$sample_id) || !nzchar(x$sample_id)) {
    .ichor_abort("sample_id must be a non-empty scalar.", "ichorviz_validation_error")
  }
  if (!x$genome_build %in% c("hg19", "hg38")) {
    .ichor_abort("genome_build must be hg19 or hg38.", "ichorviz_validation_error")
  }
  bins <- x$bins
  if (!nrow(bins)) .ichor_abort("The sample has no bins.", "ichorviz_validation_error")
  bad_chr <- setdiff(unique(bins$chr), .chr_levels)
  if (length(bad_chr)) {
    .ichor_abort(sprintf("Unsupported chromosome%s: %s",
                         if (length(bad_chr) > 1) "s" else "",
                         paste(bad_chr, collapse = ", ")),
                 "ichorviz_validation_error")
  }
  if (any(!is.finite(bins$start) | !is.finite(bins$end) |
          bins$start < 1 | bins$end < bins$start)) {
    .ichor_abort("Bins contain invalid coordinates.", "ichorviz_validation_error")
  }
  key <- paste(bins$chr, bins$start, bins$end, sep = ":")
  if (anyDuplicated(key)) .ichor_abort("Bins contain duplicate intervals.", "ichorviz_validation_error")
  by_chr <- split(bins, bins$chr)
  overlaps <- vapply(by_chr, function(d) {
    d <- d[order(d$start, d$end), , drop = FALSE]
    nrow(d) > 1L && any(d$start[-1] <= d$end[-nrow(d)])
  }, logical(1))
  if (any(overlaps)) {
    .ichor_abort(sprintf("Bins overlap on chromosome%s %s.",
                         if (sum(overlaps) > 1) "s" else "",
                         paste(names(overlaps)[overlaps], collapse = ", ")),
                 "ichorviz_validation_error")
  }
  if (!is.null(x$params)) {
    tf <- x$params$tumor_fraction[1]
    if (!is.na(tf) && (tf < 0 || tf > 1)) {
      .ichor_abort("Tumor fraction must be expressed as a fraction in [0, 1].",
                   "ichorviz_validation_error")
    }
    ploidy <- x$params$ploidy[1]
    if (!is.na(ploidy) && ploidy <= 0) {
      .ichor_abort("Ploidy must be positive.", "ichorviz_validation_error")
    }
  }
  invisible(x)
}

#' @export
print.ichor_sample <- function(x, ...) {
  tf <- if (is.null(x$params)) NA_real_ else x$params$tumor_fraction[1]
  cat(sprintf("<ichor_sample> %s\n", x$sample_id))
  cat(sprintf("  genome:   %s\n", x$genome_build))
  cat(sprintf("  bins:     %s\n", format(nrow(x$bins), big.mark = ",")))
  cat(sprintf("  segments: %s\n", if (is.null(x$segments)) "not loaded" else nrow(x$segments)))
  cat(sprintf("  TF:       %s\n", if (is.na(tf)) "not available" else sprintf("%.3f", tf)))
  invisible(x)
}
