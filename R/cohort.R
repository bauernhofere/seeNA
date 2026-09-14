.resolve_manifest_path <- function(path, root) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  expanded <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) {
    if (is.null(root)) .ichor_abort("Relative paths in a data-frame manifest require root.", "ichorviz_manifest_error")
    expanded <- file.path(root, expanded)
  }
  normalizePath(expanded, mustWork = FALSE)
}

#' Read a cohort manifest
#'
#' Reads every sample listed in a CSV/TSV manifest or data frame into a
#' validated `ichor_cohort`, keeping manifest order and extra columns as
#' sample metadata.
#'
#' @details
#' Required columns: `sample_id`, `cna_seg`; optional: `seg`, `params`. Other
#' columns become plot annotations; use `metadata_columns` to allowlist them.
#' Manifest IDs are aliases; source identities in the files are still checked.
#' Loading is atomic: any failed sample aborts the whole import with one
#' message listing every failure.
#' @param manifest CSV/TSV path or data frame.
#' @param genome_build Explicit shared build, `"hg19"` or `"hg38"`.
#' @param workers Positive integer; Unix uses fork workers, Windows sequential.
#' @param root Required for relative paths in a data-frame manifest. File
#'   manifests resolve paths against their own directory.
#' @param retain_paths Retain absolute paths in sample provenance? Default `FALSE`.
#' @param bounds Interval bounds policy passed to [read_ichor_sample()].
#' @param metadata_columns Optional annotation allowlist; `NULL` retains all
#'   non-file columns, `character()` retains only `sample_id`.
#' @return A validated `ichor_cohort` with `samples` and `metadata` in
#'   manifest order.
#' @examples
#' root <- system.file("extdata", package = "ichorViz")
#' cohort <- read_ichor_cohort(file.path(root, "example-manifest.csv"), "hg38")
#' cohort
#' cohort$metadata
#' @export
read_ichor_cohort <- function(manifest, genome_build, workers = 1L, root = NULL,
                              retain_paths = FALSE, bounds = c("window", "error", "trim"),
                              metadata_columns = NULL) {
  genome_build <- match.arg(genome_build, c("hg19", "hg38"))
  bounds <- match.arg(bounds)
  .scalar(workers, "workers", integer = TRUE, lower = 1, upper = .Machine$integer.max)
  .flag(retain_paths, "retain_paths")
  if (is.character(manifest) && length(manifest) == 1L) {
    path <- .assert_file(manifest, "manifest")
    root <- dirname(path)
    tab <- data.table::fread(file = path, data.table = FALSE, check.names = FALSE,
                            colClasses = list(character = "sample_id"))
  } else if (is.data.frame(manifest)) {
    tab <- as.data.frame(manifest, stringsAsFactors = FALSE)
  } else .ichor_abort("manifest must be a file path or data frame.", "ichorviz_manifest_error")
  if (!all(c("sample_id", "cna_seg") %in% names(tab)) || anyDuplicated(names(tab))) {
    .ichor_abort("Manifest needs unique columns including sample_id and cna_seg.", "ichorviz_manifest_error")
  }
  tab$sample_id <- as.character(tab$sample_id)
  if (!nrow(tab) || anyNA(tab$sample_id) || any(!nzchar(trimws(tab$sample_id))) || anyDuplicated(tab$sample_id)) {
    .ichor_abort("Manifest IDs must be non-empty and unique.", "ichorviz_manifest_error")
  }
  file_cols <- intersect(c("cna_seg", "seg", "params"), names(tab))
  for (nm in file_cols) tab[[nm]] <- vapply(as.character(tab[[nm]]), .resolve_manifest_path, character(1), root = root)
  available <- setdiff(names(tab), file_cols)
  if (!is.null(metadata_columns) && !all(metadata_columns %in% available)) .ichor_abort("Unknown metadata columns.")
  # NULL keeps every non-file column; character() is an explicit empty allowlist.
  keep <- if (is.null(metadata_columns)) available else metadata_columns
  metadata <- tab[unique(c("sample_id", keep))]
  load_one <- function(i) {
    tryCatch(withCallingHandlers(list(sample = read_ichor_sample(tab$cna_seg[i],
      seg = if ("seg" %in% names(tab)) tab$seg[i] else NULL,
      params = if ("params" %in% names(tab)) tab$params[i] else NULL,
      genome_build = genome_build, sample_id = tab$sample_id[i],
      retain_paths = retain_paths, bounds = bounds)),
      ichorviz_bounds_warning = function(w) invokeRestart("muffleWarning"),
      ichorviz_window_padding = function(m) invokeRestart("muffleMessage")),
      error = function(e) list(error = conditionMessage(e)))
  }
  if (workers > 1L && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(nrow(tab)), load_one, mc.cores = min(workers, nrow(tab)))
  } else {
    if (workers > 1L) warning("Windows uses sequential file loading.", call. = FALSE)
    results <- lapply(seq_len(nrow(tab)), load_one)
  }
  ok <- vapply(results, function(r) is.list(r) && inherits(r$sample, "ichor_sample"), logical(1))
  if (!all(ok)) {
    details <- vapply(which(!ok), function(i) {
      reason <- if (is.list(results[[i]])) results[[i]]$error else "Worker failure"
      sprintf("row %d (%s): %s", i, tab$sample_id[i], .default_if_null(reason, "Worker failure"))
    }, character(1))
    .ichor_abort(paste("Cohort import aborted:", paste(details, collapse = "; ")), "ichorviz_manifest_error")
  }
  samples <- lapply(results, `[[`, "sample")
  names(samples) <- tab$sample_id
  x <- structure(list(schema_version = 1L, samples = samples, metadata = metadata,
                      genome_build = genome_build), class = "ichor_cohort")
  validate_ichor_cohort(x)
  changes <- vapply(samples, function(s) if (is.null(s$coordinate_changes)) 0L else nrow(s$coordinate_changes), integer(1))
  .bounds_notice(bounds, sum(changes), sum(changes > 0))
  x
}

#' Validate a cohort
#'
#' Checks that every sample validates, identifiers are unique, metadata rows
#' match sample order and all samples share the cohort genome build.
#' @param x An `ichor_cohort`.
#' @return `x`, invisibly.
#' @examples
#' root <- system.file("extdata", package = "ichorViz")
#' cohort <- read_ichor_cohort(file.path(root, "example-manifest.csv"), "hg38")
#' validate_ichor_cohort(cohort)
#' @export
validate_ichor_cohort <- function(x) {
  if (!inherits(x, "ichor_cohort") || !identical(x$schema_version, 1L) ||
      !is.list(x$samples) || !length(x$samples) || !is.data.frame(x$metadata)) .ichor_abort("Invalid cohort object.")
  lapply(x$samples, validate_ichor_sample)
  ids <- vapply(x$samples, `[[`, character(1), "sample_id")
  if (anyDuplicated(ids) || !identical(names(x$samples), unname(ids)) ||
      !identical(x$metadata$sample_id, unname(ids))) .ichor_abort("Cohort sample and metadata order disagree.")
  if (!all(vapply(x$samples, function(s) identical(s$genome_build, x$genome_build), logical(1)))) .ichor_abort("Cohort genome builds disagree.")
  invisible(x)
}

#' @export
print.ichor_cohort <- function(x, ...) {
  cat(sprintf("<ichor_cohort>\n  samples: %d\n  genome: %s\n", length(x$samples), x$genome_build))
  invisible(x)
}
