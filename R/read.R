.match_columns <- function(nms, specification, required) {
  keys <- .normalize_key(nms)
  found <- stats::setNames(rep(NA_integer_, length(specification)), names(specification))
  for (target in names(specification)) {
    aliases <- unique(.normalize_key(c(target, specification[[target]])))
    candidates <- which(keys %in% aliases)
    if (length(candidates) > 1L) {
      .ichor_abort(paste("Ambiguous columns for", target), "ichorviz_schema_error")
    }
    if (length(candidates)) found[target] <- candidates
  }
  if (anyDuplicated(found[!is.na(found)])) .ichor_abort("Ambiguous semantic columns.", "ichorviz_schema_error")
  missing <- required[is.na(found[required])]
  if (length(missing)) {
    .ichor_abort(paste("Missing required columns:", paste(missing, collapse = ", ")),
                 "ichorviz_schema_error")
  }
  found
}

.read_table <- function(path) {
  # No heuristic shell-command dispatch in fread; every input is a literal file.
  out <- withCallingHandlers(
    data.table::fread(file = path, data.table = FALSE, check.names = FALSE,
                     colClasses = "character", na.strings = c("NA", "NaN", "")),
    warning = function(w) .ichor_abort("Malformed or empty tabular input.", "ichorviz_schema_error")
  )
  if (!nrow(out) || anyDuplicated(names(out))) {
    .ichor_abort("Empty table or duplicate column names.", "ichorviz_schema_error")
  }
  out
}

.normalize_table <- function(raw, spec, required) {
  idx <- .match_columns(names(raw), spec, required)
  keep <- idx[!is.na(idx)]
  out <- raw[unname(keep)]
  names(out) <- names(keep)
  out$chr <- .normalize_chr(out$chr)
  for (nm in intersect(c("start", "end", "logR", "median", "bins", "copy_number",
                         "corrected_copy_number", "logR_copy_number"), names(out))) {
    # ichorCNA's diagnostic inverse transform can be infinite at zero TF.
    out[[nm]] <- .as_number(out[[nm]], nm, allow_infinite = nm == "logR_copy_number")
  }
  for (nm in intersect(c("event", "corrected_call"), names(out))) out[[nm]] <- .normalize_call(out[[nm]])
  if ("subclone_status" %in% names(out)) out$subclone_status <- .as_flag(out$subclone_status)
  out <- out[order(.chr_rank(out$chr), out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Read an ichorCNA bin-level file
#'
#' Accepts one sample per file. Mixed sample prefixes and duplicate semantic
#' columns are errors. Missing scientific values remain NA; malformed values
#' are errors. Coordinates are 1-based closed intervals.
#' @param path Path to a `.cna.seg` file.
#' @param sample_id Optional output alias, not a selector for multi-sample files.
#' @return Normalized data frame with `sample_id` and `source_id` attributes.
#' @export
read_ichor_cna <- function(path, sample_id = NULL) {
  path <- .assert_file(path, "cna_seg")
  raw <- .read_table(path)
  spec <- list(chr = c("chromosome", "chrom"), start = "chrom_start",
               end = c("stop", "chrom_end"), corrected_copy_number = "corrected_cn",
               corrected_call = "corrected_event", copy_number = "cn", event = "call",
               logR = c("log_ratio", "log2_ratio"), subclone_status = "subclone",
               logR_copy_number = "logr_cn")
  aliases <- unique(.normalize_key(c(names(spec), unlist(spec))))
  fields <- names(raw)
  prefixes <- character()
  for (i in seq_along(fields)) {
    if (.normalize_key(fields[i]) %in% aliases) next
    # Strip only a literal dot-separated sample prefix, never an arbitrary suffix.
    dots <- gregexpr(".", fields[i], fixed = TRUE)[[1]]
    dots <- dots[dots > 0]
    hits <- dots[vapply(dots, function(pos) {
      .normalize_key(substring(fields[i], pos + 1)) %in% aliases
    }, logical(1))]
    if (length(hits)) {
      pos <- min(hits)
      prefixes <- c(prefixes, substr(fields[i], 1, pos - 1))
      fields[i] <- substring(fields[i], pos + 1)
    }
  }
  if (length(unique(prefixes)) > 1L) .ichor_abort("Multiple sample prefixes in .cna.seg.", "ichorviz_schema_error")
  names(raw) <- fields
  id <- if (length(prefixes)) prefixes[1] else sub("\\.cna\\.seg$", "", basename(path))
  out <- .normalize_table(raw, spec, c("chr", "start", "end", "logR"))
  attr(out, "source_id") <- id
  attr(out, "sample_id") <- sample_id %||% id
  out
}

#' Read an ichorCNA segment file
#'
#' Supports v0.4 `.seg` and its separate `.seg.txt` export. The latter uses
#' `chrom`, `num.mark`, and `seg.median.logR` headers.
#' @param path Path to a segment file.
#' @return A data frame with normalized segment columns.
#' @export
read_ichor_segments <- function(path) {
  raw <- .read_table(.assert_file(path, "seg"))
  spec <- list(sample_id = c("sample", "id"), chr = c("chromosome", "chrom"),
               start = c("chrom_start", "loc_start"), end = c("stop", "chrom_end", "loc_end"),
               event = "call", copy_number = "cn", bins = c("n_bins", "num_bins", "num_mark"),
               median = c("segment_median", "seg_mean", "seg_median_logR"),
               subclone_status = "subclone")
  out <- .normalize_table(raw, spec, c("chr", "start", "end", "median"))
  if ("sample_id" %in% names(out) && (anyNA(out$sample_id) || length(unique(out$sample_id)) != 1L)) {
    .ichor_abort("Segment file must contain exactly one source sample ID.", "ichorviz_schema_error")
  }
  attr(out, "source_id") <- if ("sample_id" %in% names(out)) out$sample_id[1] else
    sub("\\.seg(\\.txt)?$", "", basename(path))
  out
}

#' Read an ichorCNA parameter file
#'
#' Accepts single-sample tabular and/or key-value sections. Repeated scientific
#' fields must agree numerically. Empty files and multiple samples are errors.
#' @param path Path to a `.params.txt` file.
#' @return One-row data frame: `sample_id`, `tumor_fraction`, `ploidy`, `gender`.
#' @export
read_ichor_params <- function(path) {
  path <- .assert_file(path, "params")
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  values <- list()
  add <- function(key, val) {
    if (key %in% c("tumor_fraction", "ploidy")) val <- .as_number(val, key) else val <- .text_missing(val)
    if (key %in% names(values) && !isTRUE(all.equal(values[[key]], val, tolerance = 1e-8))) {
      .ichor_abort(paste("Conflicting parameter values for", key), "ichorviz_schema_error")
    }
    values[[key]] <<- val
  }
  header <- grep("^Sample\\t", lines, ignore.case = TRUE)
  if (length(header)) {
    if (length(header) != 1 || header != 1 || length(lines) < 2) {
      .ichor_abort("Malformed parameter header.", "ichorviz_schema_error")
    }
    keys <- .normalize_key(strsplit(lines[1], "\t", fixed = TRUE)[[1]])
    vals <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
    if (length(keys) != length(vals) || anyDuplicated(keys)) .ichor_abort("Malformed parameter row.", "ichorviz_schema_error")
    for (i in seq_along(keys)) if (keys[i] %in% c("sample", "tumor_fraction", "ploidy", "gender")) add(keys[i], vals[i])
    lines <- lines[-(1:2)]
  }
  # v0.4 appends a candidate-initialization table after the chosen solution.
  # It is diagnostic, not another selected sample or a parameter override.
  diagnostics <- grep("^init\\tn_est\\tphi_est\\t", lines)
  if (length(diagnostics)) lines <- lines[seq_len(diagnostics[1] - 1L)]
  for (line in lines) {
    if (!grepl(":", line, fixed = TRUE)) {
      if (grepl("\t", line, fixed = TRUE)) .ichor_abort("Multiple sample rows in parameter file.", "ichorviz_schema_error")
      add("sample", line)
      next
    }
    key <- .normalize_key(sub(":.*$", "", line))
    if (key %in% c("sample", "tumor_fraction", "ploidy", "gender")) add(key, trimws(sub("^[^:]*:", "", line)))
  }
  if (!any(c("tumor_fraction", "ploidy") %in% names(values))) {
    .ichor_abort("No recognized scientific parameters.", "ichorviz_schema_error")
  }
  data.frame(sample_id = values$sample %||% sub("\\.params\\.txt$", "", basename(path)),
             tumor_fraction = values$tumor_fraction %||% NA_real_,
             ploidy = values$ploidy %||% NA_real_, gender = values$gender %||% NA_character_,
             stringsAsFactors = FALSE)
}
