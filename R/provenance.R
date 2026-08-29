#' Summarize input provenance
#'
#' @param x An `ichor_sample` or `ichor_cohort`.
#' @return A data frame of source files, sizes, and MD5 checksums.
#' @export
ichor_provenance <- function(x) {
  samples <- if (inherits(x, "ichor_sample")) list(x) else if (inherits(x, "ichor_cohort")) x$samples else
    .ichor_abort("x must be an ichor_sample or ichor_cohort.")
  rows <- lapply(samples, function(s) {
    paths <- unlist(s$provenance, use.names = TRUE)
    paths <- paths[!is.na(paths) & nzchar(paths)]
    info <- file.info(paths)
    data.frame(sample_id = s$sample_id, type = names(paths), path = unname(paths),
               bytes = info$size, md5 = unname(tools::md5sum(paths)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
