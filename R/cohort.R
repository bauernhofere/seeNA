# Cohorts and aligned matrices --------------------------------------------

.resolve_manifest_path <- function(path, root) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  expanded <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) expanded <- file.path(root, expanded)
  normalizePath(expanded, mustWork = FALSE)
}

#' Read a cohort manifest
#'
#' The manifest must contain `sample_id` and `cna_seg`. Optional `seg` and
#' `params` columns identify companion files; all other columns are retained as
#' sample metadata. Relative paths are resolved from the manifest's directory.
#'
#' @param manifest Path to a CSV/TSV manifest or a data frame.
#' @param genome_build Genome build shared by the cohort.
#' @param workers Number of file-reading processes. Values above one use
#'   `parallel::mclapply()` on Unix-like systems.
#' @return An `ichor_cohort` object.
#' @export
read_ichor_cohort <- function(manifest, genome_build, workers = 1L) {
  genome_build <- match.arg(genome_build, c("hg19", "hg38"))
  if (is.character(manifest) && length(manifest) == 1L) {
    path <- .assert_file(manifest, "manifest")
    root <- dirname(path)
    tab <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else if (is.data.frame(manifest)) {
    tab <- as.data.frame(manifest, stringsAsFactors = FALSE)
    root <- getwd()
  } else {
    .ichor_abort("manifest must be a path or data frame.", "ichorviz_manifest_error")
  }
  required <- c("sample_id", "cna_seg")
  missing <- setdiff(required, names(tab))
  if (length(missing)) .ichor_abort(sprintf("Manifest is missing: %s", paste(missing, collapse = ", ")),
                                     "ichorviz_manifest_error")
  if (!nrow(tab) || anyNA(tab$sample_id) || any(!nzchar(tab$sample_id)) || anyDuplicated(tab$sample_id)) {
    .ichor_abort("Manifest sample_id values must be non-empty and unique.", "ichorviz_manifest_error")
  }
  for (nm in intersect(c("cna_seg", "seg", "params"), names(tab))) {
    tab[[nm]] <- vapply(as.character(tab[[nm]]), .resolve_manifest_path, character(1), root = root)
  }

  load_one <- function(i) {
    read_ichor_sample(
      cna_seg = tab$cna_seg[i],
      seg = if ("seg" %in% names(tab)) tab$seg[i] else NULL,
      params = if ("params" %in% names(tab)) tab$params[i] else NULL,
      genome_build = genome_build,
      sample_id = tab$sample_id[i]
    )
  }
  workers <- as.integer(workers)
  if (!is.finite(workers) || workers < 1L) .ichor_abort("workers must be a positive integer.")
  if (workers > 1L && .Platform$OS.type != "windows") {
    samples <- parallel::mclapply(seq_len(nrow(tab)), load_one, mc.cores = workers)
  } else {
    if (workers > 1L) warning("Parallel manifest loading currently falls back to sequential on Windows.", call. = FALSE)
    samples <- lapply(seq_len(nrow(tab)), load_one)
  }
  names(samples) <- tab$sample_id
  metadata <- tab[setdiff(names(tab), c("cna_seg", "seg", "params"))]
  structure(list(samples = samples, metadata = metadata, genome_build = genome_build,
                 manifest = if (exists("path")) path else NULL), class = "ichor_cohort")
}

#' @export
print.ichor_cohort <- function(x, ...) {
  cat("<ichor_cohort>\n")
  cat(sprintf("  samples: %s\n", length(x$samples)))
  cat(sprintf("  genome:  %s\n", x$genome_build))
  invisible(x)
}

.rebin_sample <- function(sample, bin_size, value) {
  d <- sample$bins
  if (value == "call") {
    call_col <- if ("corrected_call" %in% names(d)) "corrected_call" else "event"
    if (is.null(call_col) || !call_col %in% names(d)) .ichor_abort("No call column is available.")
    d$value <- .call_score(d[[call_col]])
  } else {
    if (!value %in% names(d)) .ichor_abort(sprintf("%s is unavailable for %s.", value, sample$sample_id))
    d$value <- .as_number(d[[value]])
  }
  d <- d[!is.na(d$value), c("chr", "start", "end", "value"), drop = FALSE]
  if (!nrow(d)) return(data.frame(chr = character(), bin_index = numeric(), value = numeric()))

  first <- floor((d$start - 1) / bin_size)
  last <- floor((d$end - 1) / bin_size)
  counts <- last - first + 1
  source_row <- rep.int(seq_len(nrow(d)), counts)
  bin_index <- if (all(counts == 1)) first else
    unlist(Map(seq.int, first, last), use.names = FALSE)
  target_start <- bin_index * bin_size + 1
  target_end <- (bin_index + 1) * bin_size
  overlap <- pmax(0, pmin(d$end[source_row], target_end) -
                    pmax(d$start[source_row], target_start) + 1)
  chr <- d$chr[source_row]
  weighted <- d$value[source_row] * overlap
  key <- paste(chr, bin_index, sep = ":")

  if (!anyDuplicated(key)) {
    return(data.frame(chr = chr, bin_index = bin_index, value = weighted / overlap))
  }
  unique_key <- unique(key)
  group <- match(key, unique_key)
  first_row <- match(unique_key, key)
  numerator <- as.numeric(rowsum(weighted, group, reorder = FALSE))
  denominator <- as.numeric(rowsum(overlap, group, reorder = FALSE))
  data.frame(chr = chr[first_row], bin_index = bin_index[first_row],
             value = numerator / denominator)
}

#' Build an aligned cohort bin matrix
#'
#' Source bins are assigned to fixed genomic bins using overlap-weighted values.
#'
#' @param cohort An `ichor_cohort`.
#' @param bin_size Target bin width in base pairs.
#' @param value One of `"logR"`, `"corrected_copy_number"`, `"copy_number"`,
#'   or `"call"`. Calls are encoded as deep loss `-2`, loss `-1`, neutral `0`,
#'   and gain `1`.
#' @return An `ichor_matrix` containing `values`, genomic `bins`, and sample
#'   metadata.
#' @export
ichor_matrix <- function(cohort, bin_size = 1e6,
                         value = c("logR", "corrected_copy_number", "copy_number", "call")) {
  if (!inherits(cohort, "ichor_cohort")) .ichor_abort("cohort must be an ichor_cohort.")
  value <- match.arg(value)
  bin_size <- as.numeric(bin_size)
  if (!is.finite(bin_size) || bin_size <= 0) .ichor_abort("bin_size must be positive.")

  sizes <- .chromosome_sizes(cohort$genome_build)
  grid <- do.call(rbind, lapply(seq_len(nrow(sizes)), function(i) {
    idx <- 0:(ceiling(sizes$length[i] / bin_size) - 1)
    data.frame(chr = sizes$chr[i], bin_index = idx,
               start = idx * bin_size + 1,
               end = pmin((idx + 1) * bin_size, sizes$length[i]))
  }))
  grid$key <- paste(grid$chr, grid$bin_index, sep = ":")
  grid$label <- sprintf("chr%s:%s-%s", grid$chr, grid$start, grid$end)

  mat <- matrix(NA_real_, nrow = length(cohort$samples), ncol = nrow(grid),
                dimnames = list(names(cohort$samples), grid$label))
  for (i in seq_along(cohort$samples)) {
    d <- .rebin_sample(cohort$samples[[i]], bin_size, value)
    key <- paste(d$chr, d$bin_index, sep = ":")
    j <- match(key, grid$key)
    mat[i, j[!is.na(j)]] <- d$value[!is.na(j)]
  }
  grid$key <- NULL
  structure(list(values = mat, bins = grid, samples = cohort$metadata,
                 genome_build = cohort$genome_build, bin_size = bin_size, value = value),
            class = "ichor_matrix")
}

#' @export
print.ichor_matrix <- function(x, ...) {
  cat("<ichor_matrix>\n")
  cat(sprintf("  dimensions: %s samples x %s bins\n", nrow(x$values), ncol(x$values)))
  cat(sprintf("  value:      %s\n", x$value))
  cat(sprintf("  bin size:   %s bp\n", format(x$bin_size, big.mark = ",", scientific = FALSE)))
  invisible(x)
}
