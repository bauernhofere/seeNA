.check_window_padding <- function(bins, segments, build) {
  sizes <- .chromosome_sizes(build)
  limit <- sizes$length[match(bins$chr, sizes$chr)]
  over <- which(bins$end > limit)
  seg_limit <- sizes$length[match(segments$chr, sizes$chr)]
  seg_over <- which(segments$end > seg_limit)
  if (!length(over) && !length(seg_over)) return(invisible(NULL))
  fail <- function() .ichor_abort(paste(
    "Out-of-build intervals do not match supported fixed-window terminal padding.",
    "Verify the genome build and source bin grid; use bounds = 'error' for strict bounds",
    "or explicit bounds = 'trim' only after inspecting the original intervals."), "ichorviz_bounds_error")

  # Evidence comes from observed interior bins, never the amount of overhang alone.
  interior <- which(is.finite(bins$start) & is.finite(bins$end) & bins$start >= 1 & bins$end < limit)
  widths <- unique(bins$end[interior] - bins$start[interior] + 1)
  if (length(interior) < 2L || length(widths) != 1L || !is.finite(widths) || widths < 1 || widths != floor(widths)) fail()
  w <- widths
  grid_ok <- is.finite(bins$start) & is.finite(bins$end) & !is.na(limit) &
    bins$start >= 1 & bins$start <= limit & (bins$start - 1) %% w == 0 &
    (bins$end == bins$start + w - 1 | bins$end == pmin(bins$start + w - 1, limit))
  if (!all(grid_ok)) fail()
  if (length(over) && any(bins$end[over] != ceiling(limit[over] / w) * w |
                          bins$end[over] - limit[over] >= w)) fail()
  if (length(seg_over)) {
    observed <- paste(bins$chr[over], bins$end[over], sep = ":")
    if (any(!is.finite(segments$start[seg_over]) | segments$start[seg_over] > seg_limit[seg_over]) ||
        !all(paste(segments$chr[seg_over], segments$end[seg_over], sep = ":") %in% observed)) fail()
  }
  invisible(w)
}

.bounds_notice <- function(bounds, n_intervals, n_samples = 1L) {
  if (!n_intervals) return(invisible(NULL))
  text <- sprintf("Trimmed %d terminal intervals across %d sample(s) under bounds = '%s'; see coordinate_changes. Genome build remains caller-specified.",
                  n_intervals, n_samples, bounds)
  if (bounds == "trim") {
    warning(structure(list(message = text, call = NULL),
                      class = c("ichorviz_bounds_warning", "warning", "condition")))
  } else {
    message(structure(list(message = paste0(text, "\n"), call = NULL),
                      class = c("ichorviz_window_padding", "message", "condition")))
  }
  invisible(NULL)
}
