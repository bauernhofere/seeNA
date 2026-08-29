`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x))) y else x

.ichor_abort <- function(message, class = "ichorviz_error") {
  cond <- structure(list(message = message, call = NULL),
                    class = c(class, "error", "condition"))
  stop(cond)
}

.assert_file <- function(path, label) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    .ichor_abort(sprintf("%s does not exist: %s", label, path), "ichorviz_file_error")
  }
  normalizePath(path, mustWork = TRUE)
}

.as_number <- function(x) suppressWarnings(as.numeric(x))

.as_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

.normalize_key <- function(x) {
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", x)
}

.normalize_chr <- function(x) {
  x <- toupper(sub("^CHR", "", trimws(as.character(x)), ignore.case = TRUE))
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x
}

.chr_levels <- c(as.character(1:22), "X", "Y")

.chr_rank <- function(x) match(.normalize_chr(x), .chr_levels)

.call_state <- function(x) {
  x <- toupper(as.character(x))
  out <- rep("Neutral", length(x))
  out[x %in% c("GAIN", "AMP") | grepl("^HLAMP", x)] <- "Gain"
  out[x == "HETD"] <- "Loss"
  out[x == "HOMD"] <- "Deep loss"
  factor(out, levels = c("Deep loss", "Loss", "Neutral", "Gain"))
}

.call_score <- function(x) {
  s <- as.character(.call_state(x))
  unname(c("Deep loss" = -2, "Loss" = -1, "Neutral" = 0, "Gain" = 1)[s])
}
