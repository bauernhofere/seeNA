# Individual and comparison plots -----------------------------------------

.bin_state <- function(d, call_column) {
  if (!call_column %in% names(d)) {
    available <- intersect(c("corrected_call", "event"), names(d))
    if (!length(available)) return(factor(rep("Neutral", nrow(d)), levels = names(ichor_state_colors())))
    call_column <- available[1]
  }
  .call_state(d[[call_column]])
}

.segment_state <- function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  if ("event" %in% names(d)) .call_state(d$event) else
    factor(rep("Neutral", nrow(d)), levels = names(ichor_state_colors()))
}

#' Plot one genome-wide ichorCNA profile
#'
#' @param x An `ichor_sample`.
#' @param call_column Column used to color copy-number states.
#' @param colors Named state color vector.
#' @param point_size Bin point size.
#' @return A `ggplot` object.
#' @export
plot_ichor_profile <- function(x, call_column = "corrected_call",
                               colors = ichor_state_colors(), point_size = 0.35) {
  validate_ichor_sample(x)
  layout <- .genome_layout(x$genome_build, x$bins$chr)
  bins <- .add_genome_coordinates(x$bins, layout)
  bins$state <- .bin_state(bins, call_column)
  seg <- x$segments
  if (!is.null(seg)) {
    seg <- .add_genome_coordinates(seg, layout, segments = TRUE)
    seg$state <- .segment_state(seg)
  }

  p <- ggplot2::ggplot(bins, ggplot2::aes(x = x, y = logR)) +
    ggplot2::geom_vline(xintercept = utils::head(layout$boundary, -1), color = "grey90", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(color = state), size = point_size, alpha = 0.85,
                        na.rm = TRUE) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = layout$mid, labels = layout$chr,
                                expand = ggplot2::expansion(mult = c(0.002, 0.002))) +
    ggplot2::labs(y = expression(log[2]~ratio), color = "Copy-number state",
                  title = x$sample_id) +
    .theme_ichor()
  if (!is.null(seg) && nrow(seg)) {
    p <- p + ggplot2::geom_segment(
      data = seg, ggplot2::aes(x = x, xend = xend, y = median, yend = median, color = state),
      inherit.aes = FALSE, linewidth = 0.55, lineend = "round", na.rm = TRUE
    )
  }
  p
}

.as_sample_list <- function(samples) {
  if (inherits(samples, "ichor_sample")) samples <- list(samples)
  if (!is.list(samples) || !length(samples) ||
      !all(vapply(samples, inherits, logical(1), what = "ichor_sample"))) {
    .ichor_abort("samples must be an ichor_sample or a list of ichor_sample objects.")
  }
  ids <- vapply(samples, `[[`, character(1), "sample_id")
  if (anyDuplicated(ids)) .ichor_abort("Sample identifiers must be unique.")
  names(samples) <- ids
  samples
}

.comparison_data <- function(samples, region = NULL) {
  builds <- unique(vapply(samples, `[[`, character(1), "genome_build"))
  if (length(builds) != 1L) .ichor_abort("All samples must use the same genome build.")
  if (is.null(region)) {
    chromosomes <- unique(unlist(lapply(samples, function(s) s$bins$chr)))
    layout <- .genome_layout(builds, chromosomes)
  } else {
    region <- parse_ichor_region(region, builds)
    layout <- NULL
  }

  bins <- Map(function(s, id) {
    d <- s$bins
    if (is.null(region)) {
      d <- .add_genome_coordinates(d, layout)
    } else {
      d <- .subset_interval(d, region)
      d$x <- (d$start + d$end) / 2 / 1e6
    }
    d$sample <- rep(id, nrow(d))
    d
  }, samples, names(samples))
  segs <- Map(function(s, id) {
    d <- s$segments
    if (is.null(d)) return(NULL)
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
#' Overlays aligned bins and segments using one color per sample. This supports
#' paired fluids, longitudinal samples, technical replicates, and arbitrary
#' genomic regions.
#'
#' @param samples An `ichor_sample` or list of samples.
#' @param region Optional `chr:start-end` region. The default is genome-wide.
#' @param colors Optional named vector of sample colors.
#' @param point_size Bin point size.
#' @return A `ggplot` object.
#' @export
plot_ichor_compare <- function(samples, region = NULL, colors = NULL, point_size = 0.4) {
  samples <- .as_sample_list(samples)
  d <- .comparison_data(samples, region)
  if (!nrow(d$bins)) .ichor_abort("No bins overlap the requested region.", "ichorviz_region_error")
  if (is.null(colors)) colors <- stats::setNames(.default_sample_colors(length(samples)), names(samples))
  if (!all(names(samples) %in% names(colors))) .ichor_abort("colors must be named for every sample.")

  p <- ggplot2::ggplot(d$bins, ggplot2::aes(x = x, y = logR, color = sample)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
    ggplot2::geom_point(size = point_size, alpha = 0.65, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = colors[names(samples)]) +
    ggplot2::labs(y = expression(log[2]~ratio), color = NULL) +
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
                                  expand = ggplot2::expansion(mult = c(0.002, 0.002)))
  } else {
    p <- p +
      ggplot2::coord_cartesian(xlim = c(d$region$start, d$region$end) / 1e6) +
      ggplot2::labs(x = sprintf("chr%s position (Mb)", d$region$chr),
                    title = sprintf("chr%s:%s-%s", d$region$chr,
                                    format(d$region$start, big.mark = ",", scientific = FALSE),
                                    format(d$region$end, big.mark = ",", scientific = FALSE))) +
      ggplot2::theme(axis.title.x = ggplot2::element_text())
  }
  p
}

#' Plot a genomic region
#'
#' @inheritParams plot_ichor_compare
#' @return A `ggplot` object.
#' @export
plot_ichor_region <- function(samples, region, colors = NULL, point_size = 0.55) {
  plot_ichor_compare(samples, region = region, colors = colors, point_size = point_size)
}
