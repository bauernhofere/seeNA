.concordance_colors <- function() c(a_only = "#007A87", b_only = "#B34E18",
  concordant = "#35864C", discordant = "#9254A1", unavailable = "#D9D9D9",
  baseline_flag = "#D9AB38")

#' Overlay a pair with a directional call-agreement track
#'
#' Draws original bin/segment profiles above a shared-grid agreement track.
#' This is not a plot of a minus b: bar colors describe call direction, while
#' heights show selected raw or ploidy-adjusted logR. Both panels share `ylim`.
#'
#' @details
#' Default `height="both"` shows a in the left half and b in the right half of
#' each target bin when both are altered. For one-sided calls only the altered
#' sample contributes a bar. Half-bin placement is display packing, not distinct
#' sub-bin measurements. `"representative"` uses the explicit mean/largest-
#' absolute rule in [ichor_pair_concordance()], including its order-dependent tie.
#'
#' Both-neutral bins are blank. Grey marks unknown comparisons or unavailable
#' selected heights; gold strips flag unresolved X/Y neutral-reference evidence.
#' Thus unavailable evidence cannot disappear as apparently neutral white space.
#' An exactly zero altered height is shown as a colored marker, without nudging it.
#' Gold does not invalidate a source call under the default sex-chromosome policy.
#' Fixed view limits clip only the drawing and warn on out-of-view heights.
#' @inheritParams ichor_pair_concordance
#' @param region Optional chromosome or `chr:start-end`. Target bins are aggregated
#'   before regional cropping; upper points retain original source midpoints.
#' @param height `"both"` (default) or `"representative"` bar-height policy.
#' @param ylim Finite increasing two-element vector straddling zero. Default -2,2.
#' @param sample_labels Two distinct display names for a,b. Default A,B; pairing
#'   and fluid identities are never inferred. Use e.g. c("Plasma", "Urine").
#' @param sample_colors Two named colors keyed by a,b.
#' @param colors Named agreement palette including a_only, b_only, concordant,
#'   discordant, unavailable and baseline_flag.
#' @param point_size Upper-panel bin point size.
#' @param show_sample_id Include source aliases in the subtitle? Default TRUE.
#'   FALSE keeps generic/display names and fitted parameters. Plot data and
#'   attributes still contain source identifiers; this is not de-identification.
#' @return Faceted ggplot object. Attributes `ichor_concordance` (displayed-bin
#'   table), `ichor_transform` and `ichor_view` retain data and policy settings.
#' @examples
#' root <- system.file("extdata", package = "ichorViz")
#' a <- read_ichor_sample(file.path(root, "example-a.cna.seg"), genome_build = "hg38")
#' b <- read_ichor_sample(file.path(root, "example-b.cna.seg"), genome_build = "hg38")
#' plot_ichor_concordance(a, b, region = "chr1:1-5000000", show_sample_id = FALSE)
#' @export
plot_ichor_concordance <- function(a, b, region = NULL, bin_size = 1e6,
    call_column = c("corrected_call", "event"), min_coverage = 1,
    chromosomes = c(as.character(1:22), "X"), ploidy_adjust = FALSE,
    sex_chromosomes = c("flag", "require_neutral"), max_cells = 5e7,
    height = c("both", "representative"), ylim = c(-2, 2),
    sample_labels = c("A", "B"), sample_colors = c(a = "#007A87", b = "#B34E18"),
    colors = .concordance_colors(), point_size = 0.4, show_sample_id = TRUE) {
  height <- match.arg(height)
  call_column <- match.arg(call_column)
  sex_chromosomes <- match.arg(sex_chromosomes)
  .flag(show_sample_id, "show_sample_id")
  .scalar(point_size, "point_size", lower = 0)
  if (!is.numeric(ylim) || length(ylim) != 2L || any(!is.finite(ylim)) || ylim[1] >= 0 || ylim[2] <= 0) {
    .ichor_abort("ylim must be finite, increasing and straddle zero.")
  }
  if (!is.character(sample_labels) || length(sample_labels) != 2L || anyNA(sample_labels) ||
      any(!nzchar(trimws(sample_labels))) || anyDuplicated(sample_labels)) .ichor_abort("Provide two distinct sample_labels.")
  .check_colors(sample_colors, c("a", "b"))
  .check_colors(colors, names(.concordance_colors()))
  d <- ichor_pair_concordance(a, b, bin_size, call_column, min_coverage,
    chromosomes, ploidy_adjust, sex_chromosomes, max_cells)
  settings <- attr(d, "settings")
  layout <- ichor_genome_layout(a$genome_build, settings$chromosomes)
  r <- if (is.null(region)) NULL else parse_ichor_region(region, a$genome_build)
  if (!is.null(r)) {
    d <- d[d$chr == r$chr & d$end >= r$start & d$start <= r$end, , drop = FALSE]
    if (!nrow(d)) .ichor_abort("Region is outside the selected chromosome grid.", "ichorviz_region_error")
  }
  offset <- function(chr) if (is.null(r)) layout$offset[match(chr, layout$chr)] else rep(0, length(chr))
  d$xleft <- (offset(d$chr) + d$start - 0.5) / 1e6
  d$xright <- (offset(d$chr) + d$end + 0.5) / 1e6
  panels <- c("Profiles", "Directional call agreement")
  profile <- list(); segments <- list()
  for (i in 1:2) {
    s <- list(a, b)[[i]]
    select <- function(v) {
      if (is.null(v)) return(NULL)
      v <- v[v$chr %in% settings$chromosomes, , drop = FALSE]
      if (!is.null(r)) v <- v[v$chr == r$chr & v$end >= r$start & v$start <= r$end, , drop = FALSE]
      v
    }
    v <- select(s$bins)
    profile[[i]] <- data.frame(x = (offset(v$chr) + (v$start + v$end) / 2) / 1e6,
      value = v$logR + settings$shifts[i], sample = rep(c("a", "b")[i], nrow(v)), panel = rep(panels[1], nrow(v)))
    v <- select(s$segments)
    if (!is.null(v) && nrow(v)) {
      if (!is.null(r)) v <- .subset_interval(v, r)
      segments[[i]] <- data.frame(x = (offset(v$chr) + v$start) / 1e6,
        xend = (offset(v$chr) + v$end) / 1e6, value = v$median + settings$shifts[i],
        sample = c("a", "b")[i], panel = panels[1])
    }
  }
  profile <- do.call(rbind, profile)
  segments <- do.call(rbind, segments)
  bars <- data.frame(xleft = numeric(), xright = numeric(), value = numeric(), category = character(), panel = character())
  missing_height <- rep(FALSE, nrow(d))
  if (height == "both") {
    for (s in c("a", "b")) {
      selected <- d$concordance %in% c(paste0(s, "_only"), "concordant", "discordant")
      y <- d[[paste0("logR_", s)]]
      missing_height <- missing_height | (selected & !is.finite(y))
      i <- which(selected & is.finite(y))
      center <- (d$xleft[i] + d$xright[i]) / 2
      bars <- rbind(bars, data.frame(
        xleft = if (s == "a") d$xleft[i] else center,
        xright = if (s == "a") center else d$xright[i], value = y[i],
        category = d$concordance[i], panel = rep(panels[2], length(i))))
    }
  } else {
    selected <- !d$concordance %in% c("both_neutral", "unknown")
    missing_height <- selected & !is.finite(d$representative_logR)
    i <- which(selected & !missing_height)
    bars <- data.frame(xleft = d$xleft[i], xright = d$xright[i], value = d$representative_logR[i],
                      category = d$concordance[i], panel = rep(panels[2], length(i)))
  }
  shade <- function(mask, category, low, high) {
    i <- which(mask)
    data.frame(xleft = d$xleft[i], xright = d$xright[i], ymin = rep(low, length(i)),
      ymax = rep(high, length(i)), category = rep(category, length(i)), panel = rep(panels[2], length(i)))
  }
  unavailable <- shade(d$concordance == "unknown" | missing_height, "unavailable", ylim[1], ylim[2])
  flags <- shade(d$baseline_flag, "baseline_flag", ylim[1], ylim[1] + diff(ylim) * 0.035)
  profile$panel <- factor(profile$panel, levels = panels)
  if (!is.null(segments)) segments$panel <- factor(segments$panel, levels = panels)
  bars$panel <- factor(bars$panel, levels = panels)
  unavailable$panel <- factor(unavailable$panel, levels = panels)
  flags$panel <- factor(flags$panel, levels = panels)
  guides <- data.frame(panel = factor(panels, levels = panels))
  labels <- c(a_only = paste(sample_labels[1], "only"), b_only = paste(sample_labels[2], "only"),
    concordant = "Same direction", discordant = "Opposite directions",
    unavailable = "Call/height unavailable", baseline_flag = "Sex reference unresolved")
  summary <- paste(vapply(1:2, function(i) {
    s <- list(a, b)[[i]]
    paste(sample_labels[i], if (show_sample_id) paste0("[", s$sample_id, "]") else NULL, .tf_label(s))
  }, character(1)), collapse = " | ")
  p <- ggplot2::ggplot(guides) +
    ggplot2::geom_rect(data = unavailable, ggplot2::aes(xmin = xleft, xmax = xright,
      ymin = ymin, ymax = ymax, fill = category), alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, color = "grey65", linewidth = 0.3) +
    ggplot2::geom_point(data = profile, ggplot2::aes(x = x, y = value, color = sample),
      size = point_size, alpha = 0.7, na.rm = TRUE) +
    ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xleft, xmax = xright,
      ymin = pmin(0, value), ymax = pmax(0, value), fill = category)) +
    ggplot2::geom_rect(data = flags, ggplot2::aes(xmin = xleft, xmax = xright,
      ymin = ymin, ymax = ymax, fill = category)) +
    ggplot2::geom_point(data = bars[bars$value == 0, , drop = FALSE],
      ggplot2::aes(x = (xleft + xright) / 2, y = value, fill = category),
      shape = 22, size = 0.9, stroke = 0, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = sample_colors, breaks = c("a", "b"), labels = sample_labels, name = "Sample") +
    ggplot2::scale_fill_manual(values = colors, breaks = names(labels), labels = labels,
      limits = names(labels), drop = FALSE, name = "Call agreement / evidence") +
    ggplot2::facet_grid(rows = ggplot2::vars(panel)) +
    ggplot2::labs(title = "Paired profiles and directional call agreement",
      subtitle = paste(a$genome_build, call_column, summary, sep = " | "),
      y = if (ploidy_adjust) "Ploidy-adjusted log2 ratio" else "Raw log2 ratio",
      caption = if (height == "both") paste("Track: within each target bin,", sample_labels[1],
        "left /", sample_labels[2], "right. Both neutral hidden; not a numerical difference.") else
        "Representative height: one-sided / mean / largest absolute (ties A). Not a numerical difference.") + .theme_ichor()
  if (!is.null(segments)) p <- p + ggplot2::geom_segment(data = segments,
    ggplot2::aes(x = x, xend = xend, y = value, yend = value, color = sample), linewidth = 0.55, na.rm = TRUE)
  if (is.null(r)) {
    p <- p + ggplot2::geom_vline(xintercept = utils::head(layout$boundary, -1) / 1e6,
      color = "grey85", linewidth = 0.2) +
      ggplot2::scale_x_continuous(breaks = layout$mid / 1e6, labels = layout$chr) +
      ggplot2::coord_cartesian(xlim = c(0, max(layout$boundary)) / 1e6, ylim = ylim)
  } else p <- p + ggplot2::coord_cartesian(xlim = c(r$start, r$end) / 1e6, ylim = ylim) +
    ggplot2::labs(x = sprintf("chr%s position (Mb)", r$chr)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text())
  values <- c(profile$value, segments$value, bars$value)
  clipped <- sum(is.finite(values) & (values < ylim[1] | values > ylim[2]))
  if (clipped) warning("Plotted heights exceed ylim; the view is clipped, not the input data.", call. = FALSE)
  attr(p, "ichor_concordance") <- d
  attr(p, "ichor_transform") <- settings
  attr(p, "ichor_view") <- list(region = r, ylim = ylim, height = height, clipped_heights = clipped)
  p
}
