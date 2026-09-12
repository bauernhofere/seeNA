`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x))) y else x

.ichor_abort <- function(message, class = "ichorviz_error") {
  stop(structure(list(message = message, call = NULL),
                 class = unique(c(class, "ichorviz_error", "error", "condition"))))
}

.assert_file <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || !file.exists(path) || dir.exists(path)) {
    .ichor_abort(sprintf("%s must identify an existing file.", label), "ichorviz_file_error")
  }
  normalizePath(path, mustWork = TRUE)
}

.scalar <- function(x, label, integer = FALSE, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    .ichor_abort(sprintf("%s must be a finite %s in [%s, %s].", label,
                         if (integer) "integer" else "number", lower, upper))
  }
  invisible(x)
}

.flag <- function(x, label) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) .ichor_abort(paste(label, "must be TRUE or FALSE."))
}

.text_missing <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | toupper(x) %in% c("", "NA", "NAN")] <- NA_character_
  x
}

.as_number <- function(x, label = "numeric field", allow_infinite = FALSE) {
  text <- .text_missing(x)
  value <- suppressWarnings(as.numeric(text))
  if (any(!is.na(text) & (is.na(value) | (!allow_infinite & !is.finite(value))))) {
    .ichor_abort(paste("Malformed or infinite", label), "ichorviz_schema_error")
  }
  value
}

.as_flag <- function(x) {
  text <- tolower(.text_missing(x))
  yes <- c("1", "true", "t", "yes", "y")
  no <- c("0", "false", "f", "no", "n")
  if (any(!is.na(text) & !text %in% c(yes, no))) {
    .ichor_abort("Unrecognized subclone_status token.", "ichorviz_schema_error")
  }
  out <- text %in% yes
  out[is.na(text)] <- NA
  out
}

.normalize_key <- function(x) {
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", x)
}

.normalize_chr <- function(x) {
  x <- toupper(sub("^CHR", "", trimws(as.character(x)), ignore.case = TRUE))
  x[which(x == "23")] <- "X"
  x[which(x == "24")] <- "Y"
  x
}

.chr_levels <- c(as.character(1:22), "X", "Y")
.chr_rank <- function(x) match(.normalize_chr(x), .chr_levels)

.normalize_call <- function(x) {
  x <- toupper(.text_missing(x))
  known <- x %in% c("HOMD", "HETD", "NEUT", "GAIN", "AMP", "HLAMP") |
    grepl("^HLAMP[0-9]+$", x)
  if (any(!is.na(x) & !known)) {
    .ichor_abort("Unknown copy-number call; missing calls must be NA, not neutral.",
                 "ichorviz_schema_error")
  }
  x
}

#' Map ichorCNA calls to display states
#' @param x Character vector of raw event or corrected-call tokens.
#' @return Factor with levels matching [ichor_state_colors()]. Missing values
#'   remain NA; unsupported calls raise a schema error. Does not infer new calls.
#' @export
ichor_call_state <- function(x) {
  x <- .normalize_call(x)
  out <- rep(NA_character_, length(x))
  out[which(x == "NEUT")] <- "Neutral"
  out[which(x %in% c("GAIN", "AMP") | grepl("^HLAMP", x))] <- "Gain"
  out[which(x == "HETD")] <- "Loss"
  out[which(x == "HOMD")] <- "Deep loss"
  factor(out, levels = c("Deep loss", "Loss", "Neutral", "Gain"))
}

.call_score <- function(x) {
  unname(c("Deep loss" = -2, "Loss" = -1, "Neutral" = 0, "Gain" = 1)[as.character(ichor_call_state(x))])
}
