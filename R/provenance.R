.fingerprint <- function(paths) {
  info <- file.info(paths)
  hashes <- unname(tools::md5sum(paths))
  if (anyNA(hashes) || anyNA(info$size)) .ichor_abort("Cannot fingerprint input files.", "ichorviz_file_error")
  data.frame(type = names(paths), bytes = info$size, md5 = hashes, stringsAsFactors = FALSE)
}

#' Retrieve immutable import provenance
#'
#' Returns the file fingerprints (MD5 and byte counts) and package, build and
#' schema metadata captured when the inputs were read.
#'
#' @details
#' Fingerprints are not recomputed from the current files. Paths are absent
#' unless `retain_paths = TRUE` was used at import. MD5 is a reproducibility
#' checksum, not a security or de-identification measure.
#' @param x An `ichor_sample`, `ichor_cohort`, or `ichor_matrix`.
#' @return A data frame of import-time fingerprints and metadata.
#' @examples
#' root <- system.file("extdata", package = "ichorViz")
#' a <- read_ichor_sample(file.path(root, "example-a.cna.seg"), genome_build = "hg38")
#' ichor_provenance(a)
#' @export
ichor_provenance <- function(x) {
  if (inherits(x, "ichor_sample") || inherits(x, "ichor_matrix")) return(x$provenance)
  if (!inherits(x, "ichor_cohort")) .ichor_abort("Expected a sample, cohort, or matrix.")
  out <- do.call(rbind, lapply(x$samples, ichor_provenance))
  rownames(out) <- NULL
  out
}
