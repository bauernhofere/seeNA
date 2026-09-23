# Individual and comparison plots -----------------------------------------

.bin_state <- function(d, call_column) {
  if (length(call_column) != 1 || !call_column %in% c("corrected_call", "event") || !call_column %in% names(d)) {
    .ichor_abort("Requested call_column is unavailable; choose corrected_call or event explicitly.")
  }
  ichor_call_state(d[[call_column]])
}

.check_colors <- function(colors, required) {
  if (!is.character(colors) || is.null(names(colors)) || anyDuplicated(names(colors)) ||
      !all(required %in% names(colors)) || anyNA(colors)) .ichor_abort("colors must be named for every sample or state, without duplicates.")
  tryCatch(grDevices::col2rgb(colors), error = function(e) .ichor_abort("Invalid color values."))
  invisible(colors)
}

.tf_label <- function(s) {
  if (is.null(s$params)) return("TF/ploidy unavailable")
  sprintf("TF %s; ploidy %s", if (is.na(s$params$tumor_fraction)) "NA" else sprintf("%.2f%%", 100 * s$params$tumor_fraction),
          if (is.na(s$params$ploidy)) "NA" else sprintf("%.2f", s$params$ploidy))
}

#' Plot one genome-wide ichorCNA profile
#'
#' Draws bin log ratios on a reference-length genome axis, colored by the
#' selected call column, with segment medians drawn as grey lines.
#'
#' @details
#' Segments are grey so raw segment events are not merged into the corrected
#' bin-call legend. Missing calls are grey, not neutral. Set `ploidy_adjust`
#' to reproduce the upstream plotting shift (see [ichor_adjusted_logr()]).
#' See the installed methods contract for limits.
#' @param x An `ichor_sample`.
#' @param call_column Bin column used to color states: `"corrected_call"` or
#'   `"event"`. There is no fallback between the two.
#' @param colors Named state color vector; see [ichor_state_colors()].
#' @param point_size Bin point size.
#' @param ploidy_adjust Apply [ichor_adjusted_logr()] to bins and segments?
#'   Default `FALSE` keeps raw output; `TRUE` requires TF and ploidy.
#' @param show_sample_id Show the sample identifier in the title and subtitle?
#'   Default `TRUE`. Set `FALSE` to omit identifiers from the figure.
#' @return A `ggplot` object with an `ichor_transform` metadata attribute.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' a <- read_ichor_sample(
#'   file.path(root, "example-a.cna.seg"),
#'   file.path(root, "example-a.seg"),
#'   file.path(root, "example-a.params.txt"),
#'   genome_build = "hg38"
#' )
#' plot_ichor_profile(a)
#' plot_ichor_profile(a, ploidy_adjust = TRUE, show_sample_id = FALSE)
#' @export
plot_ichor_profile <- function(x, call_column = "corrected_call",
                               colors = ichor_state_colors(), point_size = 0.35,
                               ploidy_adjust = FALSE, show_sample_id = TRUE) {
  validate_ichor_sample(x)
  .scalar(point_size, "point_size", lower = 0)
  .flag(ploidy_adjust, "ploidy_adjust")
  .flag(show_sample_id, "show_sample_id")
  .check_colors(colors, names(ichor_state_colors()))
  if (!any(is.finite(x$bins$logR))) .ichor_abort("No finite logR values to plot.")
  layout <- ichor_genome_layout(x$genome_build,
    c(as.character(1:22), "X", intersect("Y", c(x$bins$chr, x$segments$chr))))
  bins <- .add_genome_coordinates(x$bins, layout)
  bins$state <- .bin_state(bins, call_column)
  shift <- if (ploidy_adjust) .logr_shift(x) else 0
  bins$logR <- bins$logR + shift
  seg <- x$segments
  if (!is.null(seg)) {
    seg <- .add_genome_coordinates(seg, layout, segments = TRUE)
    seg$median <- seg$median + shift
  }

  p <- ggplot2::ggplot(bins, ggplot2::aes(x = x, y = logR)) +
    ggplot2::geom_vline(xintercept = utils::head(layout$boundary, -1), color = "grey90", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(color = state), size = point_size, alpha = 0.85,
                        na.rm = TRUE, show.legend = TRUE) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE, na.value = "grey65") +
    ggplot2::scale_x_continuous(breaks = layout$mid, labels = layout$chr,
                                limits = c(0, max(layout$boundary)),
                                expand = ggplot2::expansion(mult = c(0.002, 0.002))) +
    ggplot2::labs(y = if (ploidy_adjust) "Ploidy-adjusted log2 ratio" else "Raw log2 ratio", color = "Copy-number state",
                  title = if (show_sample_id) x$sample_id else NULL,
                  subtitle = paste(x$genome_build, call_column, .tf_label(x), sep = " | ")) +
    .theme_ichor()
  if (!is.null(seg) && nrow(seg)) {
    p <- p + ggplot2::geom_segment(
      data = seg, ggplot2::aes(x = x, xend = xend, y = median, yend = median),
      inherit.aes = FALSE, color = "grey35", linewidth = 0.55, lineend = "round", na.rm = TRUE
    )
  }
  attr(p, "ichor_transform") <- list(ploidy_adjust = ploidy_adjust, shift = stats::setNames(shift, x$sample_id))
  p
}

.as_sample_list <- function(samples) {
  if (inherits(samples, "ichor_sample")) samples <- list(samples)
  if (!is.list(samples) || !length(samples) ||
      !all(vapply(samples, inherits, logical(1), what = "ichor_sample"))) {
    .ichor_abort("samples must be an ichor_sample or a list of ichor_sample objects.")
  }
  lapply(samples, validate_ichor_sample)
  ids <- vapply(samples, `[[`, character(1), "sample_id")
  if (anyDuplicated(ids)) .ichor_abort("Sample identifiers must be unique.")
  names(samples) <- ids
  samples
}

.comparison_data <- function(samples, region = NULL, ploidy_adjust = FALSE) {
  builds <- unique(vapply(samples, `[[`, character(1), "genome_build"))
  if (length(builds) != 1L) .ichor_abort("All samples must use the same genome build.")
  if (is.null(region)) {
    chromosomes <- unique(unlist(lapply(samples, function(s) c(s$bins$chr, s$segments$chr))))
    layout <- ichor_genome_layout(builds, c(as.character(1:22), "X", intersect("Y", chromosomes)))
  } else {
    region <- parse_ichor_region(region, builds)
    layout <- NULL
  }

  bins <- Map(function(s, id) {
    d <- s$bins[c("chr", "start", "end", "logR")]
    if (ploidy_adjust) d$logR <- d$logR + .logr_shift(s)
    if (is.null(region)) {
      d <- .add_genome_coordinates(d, layout)
    } else {
      # Keep the original bin midpoint, not a re-centered point after cropping.
      d <- d[d$chr == region$chr & d$end >= region$start & d$start <= region$end, , drop = FALSE]
      d$x <- (d$start + d$end) / 2 / 1e6
    }
    d$sample <- rep(id, nrow(d))
    d
  }, samples, names(samples))
  segs <- Map(function(s, id) {
    d <- s$segments
    if (is.null(d)) return(NULL)
    d <- d[c("chr", "start", "end", "median")]
    if (ploidy_adjust) d$median <- d$median + .logr_shift(s)
    if (is.null(region)) {
      d <- .add_genome_coordinates(d, layout, segments = TRUE)
    } else {
      d <- .subset_interval(d, region)
      d$x <- d$start / 1e6
      d$xend <- d$end / 1e6
    }
    d$sample <- rep(id, nrow(d))
    d
  }, samples, names(samples))

  list(bins = do.call(rbind, bins), segments = do.call(rbind, Filter(Negate(is.null), segs)),
       build = builds, region = region, layout = layout)
}

#' Compare two or more ichorCNA profiles
#'
#' Overlays the original bins and segment medians of several samples in shared
#' genomic coordinates, one color per sample, genome-wide or within a region.
#'
#' @details
#' Bins are never rebinned or re-centered; only segment endpoints are clipped
#' in region views. Which samples belong together comes from your study
#' design, not from filenames. See the installed methods contract for limits.
#' @param samples An `ichor_sample` or list of samples with unique identifiers
#'   and the same genome build.
#' @param region Optional `chr:start-end` region. The default is genome-wide.
#' @param colors Optional named vector of sample colors.
#' @param point_size Bin point size.
#' @param show_sample_id Show sample identifiers in the legend and subtitle?
#'   Default `TRUE`. Set `FALSE` to hide the legend and name-bearing subtitle.
#' @inheritParams plot_ichor_profile
#' @return A `ggplot` object with an `ichor_transform` metadata attribute.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' a <- read_ichor_sample(file.path(root, "example-a.cna.seg"),
#'                        file.path(root, "example-a.seg"), genome_build = "hg38")
#' b <- read_ichor_sample(file.path(root, "example-b.cna.seg"),
#'                        file.path(root, "example-b.seg"), genome_build = "hg38")
#' plot_ichor_compare(list(a, b))
#' plot_ichor_compare(list(a, b), colors = c("example-a" = "#007A87", "example-b" = "#B34E18"))
#' @export
plot_ichor_compare <- function(samples, region = NULL, colors = NULL, point_size = 0.4,
                               ploidy_adjust = FALSE, show_sample_id = TRUE) {
  samples <- .as_sample_list(samples)
  .scalar(point_size, "point_size", lower = 0)
  .flag(ploidy_adjust, "ploidy_adjust")
  .flag(show_sample_id, "show_sample_id")
  d <- .comparison_data(samples, region, ploidy_adjust)
  if (!nrow(d$bins)) .ichor_abort("No bins overlap the requested region.", "seena_region_error")
  if (is.null(colors)) colors <- stats::setNames(.default_sample_colors(length(samples)), names(samples))
  .check_colors(colors, names(samples))
  if (!any(is.finite(d$bins$logR))) .ichor_abort("No finite logR values to plot.")
  d$bins$sample <- factor(d$bins$sample, levels = names(samples))
  if (!is.null(d$segments)) d$segments$sample <- factor(d$segments$sample, levels = names(samples))

  p <- ggplot2::ggplot(d$bins, ggplot2::aes(x = x, y = logR, color = sample)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
    ggplot2::geom_point(size = point_size, alpha = 0.65, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = colors[names(samples)], guide = if (show_sample_id) "legend" else "none") +
    ggplot2::labs(y = if (ploidy_adjust) "Ploidy-adjusted log2 ratio" else "Raw log2 ratio", color = NULL,
                  subtitle = if (show_sample_id) {
                    paste(d$build, paste(paste(names(samples), vapply(samples, .tf_label, character(1))), collapse = " | "), sep = " | ")
                  } else d$build) +
    .theme_ichor()

  if (!is.null(d$segments) && nrow(d$segments)) {
    p <- p + ggplot2::geom_segment(
      data = d$segments,
      ggplot2::aes(x = x, xend = xend, y = median, yend = median, color = sample),
      inherit.aes = FALSE, linewidth = 0.65, lineend = "round", na.rm = TRUE
    )
  }
  if (is.null(d$region)) {
    p <- p +
      ggplot2::geom_vline(xintercept = utils::head(d$layout$boundary, -1), color = "grey90", linewidth = 0.25) +
      ggplot2::scale_x_continuous(breaks = d$layout$mid, labels = d$layout$chr,
                                  limits = c(0, max(d$layout$boundary)), expand = ggplot2::expansion(mult = c(0.002, 0.002)))
  } else {
    p <- p +
      ggplot2::coord_cartesian(xlim = c(d$region$start, d$region$end) / 1e6) +
      ggplot2::labs(x = sprintf("chr%s position (Mb)", d$region$chr),
                    title = sprintf("chr%s:%s-%s", d$region$chr,
                                    format(d$region$start, big.mark = ",", scientific = FALSE),
                                    format(d$region$end, big.mark = ",", scientific = FALSE))) +
      ggplot2::theme(axis.title.x = ggplot2::element_text())
  }
  attr(p, "ichor_transform") <- list(ploidy_adjust = ploidy_adjust,
    shift = vapply(samples, function(s) if (ploidy_adjust) .logr_shift(s) else 0, numeric(1)))
  p
}

#' Plot a genomic region
#'
#' Zooms one or more samples into a chromosome or `chr:start-end` interval
#' using the original bin midpoints; a thin wrapper around [plot_ichor_compare()].
#'
#' @inheritParams plot_ichor_compare
#' @param region Required chromosome or `chr:start-end` interval.
#' @return A `ggplot` object with an `ichor_transform` metadata attribute.
#' @examples
#' root <- system.file("extdata", package = "seeNA")
#' a <- read_ichor_sample(file.path(root, "example-a.cna.seg"),
#'                        file.path(root, "example-a.seg"), genome_build = "hg38")
#' plot_ichor_region(a, "chr1:1-4000000")
#' plot_ichor_region(a, "chrX")
#' @export
plot_ichor_region <- function(samples, region, colors = NULL, point_size = 0.55,
                              ploidy_adjust = FALSE, show_sample_id = TRUE) {
  plot_ichor_compare(samples, region = region, colors = colors, point_size = point_size,
                     ploidy_adjust = ploidy_adjust, show_sample_id = show_sample_id)
}
