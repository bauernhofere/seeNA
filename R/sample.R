.validate_intervals <- function(d, build, required, label) {
  fail <- function(msg) .ichor_abort(paste(label, msg), "seena_validation_error")
  if (!is.data.frame(d) || !all(required %in% names(d)) || !nrow(d)) fail("has no rows or required columns.")
  if (!is.character(d$chr) || anyNA(d$chr) || !all(d$chr %in% .chr_levels)) fail("has unsupported chromosomes.")
  for (nm in c("start", "end")) {
    v <- d[[nm]]
    if (!is.numeric(v) || any(!is.finite(v) | v != floor(v))) fail("has non-integer or missing coordinates.")
  }
  limits <- .chromosome_sizes(build)
  if (any(d$start < 1 | d$end < d$start | d$end > limits$length[match(d$chr, limits$chr)])) {
    fail("has invalid or out-of-build coordinates. Fixed-width terminal padding may need bounds = 'window'; verify the reference build. Explicit 'trim' requires inspecting original intervals.")
  }
  if (!identical(order(.chr_rank(d$chr), d$start, d$end), seq_len(nrow(d)))) fail("must be sorted in genomic order.")
  if (nrow(d) > 1) {
    same <- d$chr[-1] == d$chr[-nrow(d)]
    if (any(same & d$start[-1] <= d$end[-nrow(d)])) fail("contains overlapping or duplicate intervals.")
  }
  for (nm in intersect(c("logR", "median", "copy_number", "corrected_copy_number", "logR_copy_number", "bins"), names(d))) {
    v <- d[[nm]]
    if (!is.numeric(v) || (nm != "logR_copy_number" && any(!is.na(v) & !is.finite(v)))) fail(paste("has invalid", nm))
    if (nm %in% c("copy_number", "corrected_copy_number", "bins") && any(v < 0, na.rm = TRUE)) fail(paste("has negative", nm))
  }
  for (nm in intersect(c("event", "corrected_call"), names(d))) .normalize_call(d[[nm]])
  if ("subclone_status" %in% names(d) && !is.logical(d$subclone_status)) fail("has non-logical subclone_status.")
}

.trim_terminal <- function(d, build, role) {
  if (is.null(d)) return(list(data = NULL, changes = NULL))
  sizes <- .chromosome_sizes(build)
  limit <- sizes$length[match(d$chr, sizes$chr)]
  hit <- which(is.finite(d$start) & is.finite(d$end) & d$start <= limit & d$end > limit)
  changes <- NULL
  if (length(hit)) {
    changes <- data.frame(role = role, row = hit, chr = d$chr[hit], start = d$start[hit],
                          original_end = d$end[hit], end = limit[hit])
    d$end[hit] <- limit[hit]
  }
  list(data = d, changes = changes)
}

#' Read one ichorCNA sample
#'
#' Reads a bin file plus optional segment and parameter files from one
#' selected ichorCNA run into a validated `ichor_sample` object.
#'
#' @details
#' Source sample IDs must agree across all supplied files before an alias is
#' applied. Intervals are validated against the stated genome build; regular
#' terminal window padding is clipped and recorded in `coordinate_changes`.
#' Import-time fingerprints are stored for reproducibility. See the installed
#' methods contract for the full coordinate and identity rules.
#' @param cna_seg Path to a `.cna.seg` file.
#' @param seg Optional segment file.
#' @param params Optional parameter file.
#' @param genome_build Explicit `"hg19"` or `"hg38"`.
#' @param sample_id Optional output alias; source identities are still checked.
#' @param retain_paths Retain absolute source paths in provenance? Default `FALSE`.
#' @param bounds Out-of-build interval policy. Default `"window"` clips only
#'   terminal padding supported by a regular source bin grid. `"error"`
#'   rejects every out-of-build interval. `"trim"` clips any straddling
#'   interval with a warning. Every change is recorded.
#' @return A validated `ichor_sample` (schema version 1) with elements
#'   `sample_id`, `genome_build`, `bins`, `segments`, `params`,
#'   `coordinate_changes` and `provenance`.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' a <- read_ichor_sample(
#'   file.path(root, "example-a.cna.seg"),
#'   file.path(root, "example-a.seg"),
#'   file.path(root, "example-a.params.txt"),
#'   genome_build = "hg38"
#' )
#' a
#' head(a$bins)
#' a$params
#' @export
read_ichor_sample <- function(cna_seg, seg = NULL, params = NULL, genome_build,
                              sample_id = NULL, retain_paths = FALSE,
                              bounds = c("window", "error", "trim")) {
  genome_build <- match.arg(genome_build, c("hg19", "hg38"))
  bounds <- match.arg(bounds)
  .flag(retain_paths, "retain_paths")
  absent <- function(p) is.null(p) || (length(p) == 1L && (is.na(p) || !nzchar(p)))
  paths <- c(cna_seg = .assert_file(cna_seg, "cna_seg"))
  if (!absent(seg)) paths <- c(paths, seg = .assert_file(seg, "seg"))
  if (!absent(params)) paths <- c(paths, params = .assert_file(params, "params"))
  before <- .fingerprint(paths)
  bins <- read_ichor_cna(paths[["cna_seg"]])
  segments <- if ("seg" %in% names(paths)) read_ichor_segments(paths[["seg"]]) else NULL
  parameters <- if ("params" %in% names(paths)) read_ichor_params(paths[["params"]]) else NULL
  source_ids <- unique(c(attr(bins, "source_id"), attr(segments, "source_id"), parameters$sample_id))
  if (length(source_ids) != 1L || is.na(source_ids) || !nzchar(source_ids)) {
    .ichor_abort("Source sample IDs disagree across bins, segments, or parameters.", "seena_identity_error")
  }
  id <- .default_if_null(sample_id, source_ids)
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(trimws(id))) .ichor_abort("sample_id must be a non-empty character scalar.")
  after <- .fingerprint(paths)
  if (!identical(before, after)) .ichor_abort("Input files changed while being read.", "seena_file_error")
  attr(bins, "source_id") <- attr(bins, "sample_id") <- NULL
  if (!is.null(segments)) {
    segments$sample_id <- NULL
    attr(segments, "source_id") <- NULL
  }
  if (!is.null(parameters)) parameters$sample_id <- id
  changes <- NULL
  if (bounds != "error") {
    if (bounds == "window") .check_window_padding(bins, segments, genome_build)
    b <- .trim_terminal(bins, genome_build, "bins")
    s <- .trim_terminal(segments, genome_build, "segments")
    bins <- b$data; segments <- s$data
    changes <- rbind(b$changes, s$changes)
  }
  before$sample_id <- id
  before$genome_build <- genome_build
  before$package_version <- as.character(utils::packageVersion("seeNA"))
  before$schema_version <- 1L
  before$coordinate_policy <- bounds
  if (retain_paths) before$path <- unname(paths)
  x <- structure(list(schema_version = 1L, sample_id = id, genome_build = genome_build,
                      bins = bins, segments = segments, params = parameters,
                      coordinate_policy = bounds, coordinate_changes = changes,
                      provenance = before), class = "ichor_sample")
  validate_ichor_sample(x)
  .bounds_notice(bounds, if (is.null(changes)) 0L else nrow(changes))
  x
}

#' Validate an ichorCNA sample
#'
#' Checks an `ichor_sample` for schema version, chromosome names, sorted
#' non-overlapping in-build intervals, numeric measurements and parameter
#' ranges.
#' @param x An `ichor_sample` object.
#' @return `x`, invisibly; invalid objects raise classed errors.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' a <- read_ichor_sample(file.path(root, "example-a.cna.seg"), genome_build = "hg38")
#' validate_ichor_sample(a)
#' @export
validate_ichor_sample <- function(x) {
  if (!inherits(x, "ichor_sample") || !identical(x$schema_version, 1L)) .ichor_abort("Expected ichor_sample schema version 1.")
  if (!is.character(x$sample_id) || length(x$sample_id) != 1L || is.na(x$sample_id) || !nzchar(x$sample_id)) {
    .ichor_abort("Invalid sample_id.", "seena_validation_error")
  }
  if (length(x$genome_build) != 1L || !x$genome_build %in% c("hg19", "hg38")) .ichor_abort("Invalid genome_build.")
  .validate_intervals(x$bins, x$genome_build, c("chr", "start", "end", "logR"), "Bins")
  if (!is.null(x$segments)) .validate_intervals(x$segments, x$genome_build, c("chr", "start", "end", "median"), "Segments")
  if (!is.null(x$params)) {
    p <- x$params
    if (!is.data.frame(p) || nrow(p) != 1L ||
        !all(c("sample_id", "tumor_fraction", "ploidy") %in% names(p)) || !identical(p$sample_id, x$sample_id)) {
      .ichor_abort("Invalid parameter table or sample identity.", "seena_validation_error")
    }
    tf <- p$tumor_fraction; ploidy <- p$ploidy
    if (!is.numeric(tf) || any(!is.na(tf) & (!is.finite(tf) | tf < 0 | tf > 1))) {
      .ichor_abort("Tumor fraction must be a fraction in [0, 1].", "seena_validation_error")
    }
    if (!is.numeric(ploidy) || any(!is.na(ploidy) & (!is.finite(ploidy) | ploidy <= 0))) {
      .ichor_abort("Ploidy must be positive and finite.", "seena_validation_error")
    }
  }
  invisible(x)
}

#' @export
print.ichor_sample <- function(x, ...) {
  tf <- if (is.null(x$params)) NA_real_ else x$params$tumor_fraction[1]
  cat(sprintf("<ichor_sample> %s\n  genome: %s\n  bins: %d\n  TF: %s\n", x$sample_id,
              x$genome_build, nrow(x$bins), if (is.na(tf)) "not available" else sprintf("%.3f", tf)))
  invisible(x)
}
