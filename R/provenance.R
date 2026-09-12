.fingerprint <- function(paths) {
  info <- file.info(paths)
  hashes <- unname(tools::md5sum(paths))
  if (anyNA(hashes) || anyNA(info$size)) .ichor_abort("Cannot fingerprint input files.", "ichorviz_file_error")
  data.frame(type = names(paths), bytes = info$size, md5 = hashes, stringsAsFactors = FALSE)
}

#' Retrieve immutable import provenance
#'
#' Returns fingerprints captured during ingestion, not hashes of current files.
#' MD5 is a reproducibility checksum, not a security or de-identification measure.
#' Source paths are absent unless explicitly retained when reading samples.
#' @param x An `ichor_sample`, `ichor_cohort`, or `ichor_matrix`.
#' @return A data frame of import-time fingerprints and package/build metadata.
#' @export
ichor_provenance <- function(x) {
  if (inherits(x, "ichor_sample") || inherits(x, "ichor_matrix")) return(x$provenance)
  if (!inherits(x, "ichor_cohort")) .ichor_abort("Expected a sample, cohort, or matrix.")
  out <- do.call(rbind, lapply(x$samples, ichor_provenance))
  rownames(out) <- NULL
  out
}
