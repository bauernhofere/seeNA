# File readers -------------------------------------------------------------

.match_columns <- function(nms, specification, required) {
  keys <- .normalize_key(nms)
  used <- rep(FALSE, length(nms))
  found <- stats::setNames(rep(NA_integer_, length(specification)), names(specification))

  for (target in names(specification)) {
    aliases <- unique(.normalize_key(c(target, specification[[target]])))
    candidates <- which(!used & keys %in% aliases)
    if (!length(candidates)) {
      suffix_hit <- vapply(keys, function(key) {
        any(vapply(aliases, function(alias) endsWith(key, paste0("_", alias)), logical(1)))
      }, logical(1))
      candidates <- which(!used & suffix_hit)
    }
    if (length(candidates)) {
      pick <- candidates[which.min(nchar(keys[candidates]))]
      found[target] <- pick
      used[pick] <- TRUE
    }
  }

  missing <- required[is.na(found[required])]
  if (length(missing)) {
    .ichor_abort(
      sprintf("Missing required column%s: %s", if (length(missing) > 1) "s" else "",
              paste(missing, collapse = ", ")),
      "ichorviz_schema_error"
    )
  }
  found
}

.infer_sample_id <- function(nms, path) {
  hit <- grep("\\.(copy\\.number|event|logR|subclone\\.status|Corrected_Copy_Number)$",
              nms, value = TRUE, ignore.case = TRUE)
  if (length(hit)) {
    inferred <- sub("\\.(copy\\.number|event|logR|subclone\\.status|Corrected_Copy_Number)$",
                    "", hit[1], ignore.case = TRUE)
    if (nzchar(inferred)) return(inferred)
  }
  sub("\\.cna\\.seg$", "", basename(path), ignore.case = TRUE)
}

#' Read an ichorCNA bin-level file
#'
#' Reads a `.cna.seg` file and normalizes both ordinary and sample-prefixed
#' column names to a stable schema.
#'
#' @param path Path to an ichorCNA `.cna.seg` file.
#' @param sample_id Optional sample identifier. By default it is inferred from
#'   prefixed columns or the filename.
#' @return A data frame with normalized bin-level columns and a `sample_id`
#'   attribute.
#' @export
read_ichor_cna <- function(path, sample_id = NULL) {
  path <- .assert_file(path, "cna_seg")
  raw <- data.table::fread(path, data.table = FALSE, check.names = FALSE,
                           na.strings = c("NA", "NaN", ""))
  if (!nrow(raw)) .ichor_abort("The .cna.seg file contains no rows.", "ichorviz_schema_error")

  spec <- list(
    chr = c("chromosome"),
    start = c("chrom_start"),
    end = c("stop", "chrom_end"),
    corrected_copy_number = c("corrected_cn"),
    corrected_call = c("corrected_event"),
    copy_number = c("cn"),
    event = c("call"),
    logR = c("log_ratio", "log2_ratio"),
    subclone_status = c("subclone"),
    logR_copy_number = c("logr_cn")
  )
  idx <- .match_columns(names(raw), spec, required = c("chr", "start", "end", "logR"))
  keep <- idx[!is.na(idx)]
  out <- raw[unname(keep)]
  names(out) <- names(keep)

  out$chr <- .normalize_chr(out$chr)
  out$start <- .as_number(out$start)
  out$end <- .as_number(out$end)
  numeric_cols <- intersect(c("logR", "copy_number", "corrected_copy_number",
                              "logR_copy_number"), names(out))
  out[numeric_cols] <- lapply(out[numeric_cols], .as_number)
  call_cols <- intersect(c("event", "corrected_call"), names(out))
  out[call_cols] <- lapply(out[call_cols], function(x) toupper(as.character(x)))
  if ("subclone_status" %in% names(out)) out$subclone_status <- .as_flag(out$subclone_status)

  out <- out[order(.chr_rank(out$chr), out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "sample_id") <- sample_id %||% .infer_sample_id(names(raw), path)
  attr(out, "source") <- path
  out
}

#' Read an ichorCNA segment file
#'
#' @param path Path to an ichorCNA `.seg` file.
#' @return A data frame with normalized segment columns.
#' @export
read_ichor_segments <- function(path) {
  path <- .assert_file(path, "seg")
  raw <- data.table::fread(path, data.table = FALSE, check.names = FALSE,
                           na.strings = c("NA", "NaN", ""))
  spec <- list(
    sample_id = c("sample", "id"),
    chr = c("chromosome"),
    start = c("chrom_start"),
    end = c("stop", "chrom_end"),
    event = c("call"),
    copy_number = c("cn"),
    bins = c("n_bins", "num_bins"),
    median = c("segment_median", "seg_mean"),
    subclone_status = c("subclone")
  )
  idx <- .match_columns(names(raw), spec,
                        required = c("chr", "start", "end", "median"))
  keep <- idx[!is.na(idx)]
  out <- raw[unname(keep)]
  names(out) <- names(keep)
  out$chr <- .normalize_chr(out$chr)
  out$start <- .as_number(out$start)
  out$end <- .as_number(out$end)
  for (nm in intersect(c("copy_number", "bins", "median"), names(out))) {
    out[[nm]] <- .as_number(out[[nm]])
  }
  if ("event" %in% names(out)) out$event <- toupper(as.character(out$event))
  if ("subclone_status" %in% names(out)) out$subclone_status <- .as_flag(out$subclone_status)
  out <- out[order(.chr_rank(out$chr), out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "source") <- path
  out
}

#' Read an ichorCNA parameter file
#'
#' Supports the tabular header and key-value sections found in ichorCNA
#' `.params.txt` files.
#'
#' @param path Path to an ichorCNA `.params.txt` file.
#' @return A one-row data frame containing `sample_id`, `tumor_fraction`,
#'   `ploidy`, and `gender` when available.
#' @export
read_ichor_params <- function(path) {
  path <- .assert_file(path, "params")
  lines <- readLines(path, warn = FALSE)
  nonempty <- which(nzchar(trimws(lines)))
  values <- list()

  header_i <- grep("^Sample\\t", lines, ignore.case = TRUE)
  if (length(header_i)) {
    i <- header_i[1]
    value_i <- nonempty[nonempty > i][1]
    if (!is.na(value_i)) {
      header <- strsplit(lines[i], "\\t", fixed = FALSE)[[1]]
      value <- strsplit(lines[value_i], "\\t", fixed = FALSE)[[1]]
      if (length(value) >= length(header)) {
        values <- as.list(stats::setNames(value[seq_along(header)], .normalize_key(header)))
      }
    }
  }

  kv <- grep("^[^:]+:\\s*", lines, value = TRUE)
  for (line in kv) {
    pieces <- strsplit(line, ":", fixed = TRUE)[[1]]
    key <- .normalize_key(pieces[1])
    val <- trimws(paste(pieces[-1], collapse = ":"))
    if (nzchar(key) && nzchar(val)) values[[key]] <- val
  }

  sample_id <- values$sample %||% values$sample_id %||%
    sub("\\.params\\.txt$", "", basename(path), ignore.case = TRUE)
  tf <- .as_number(values$tumor_fraction %||% NA_character_)
  ploidy <- .as_number(values$ploidy %||% NA_character_)
  gender <- values$gender %||% NA_character_

  data.frame(sample_id = as.character(sample_id), tumor_fraction = tf,
             ploidy = ploidy, gender = as.character(gender),
             stringsAsFactors = FALSE, check.names = FALSE)
}
