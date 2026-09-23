.annotation_palette <- function(ann) {
  lapply(ann, function(v) {
    if (is.numeric(v)) {
      observed <- v[is.finite(v)]
      limits <- if (length(observed)) range(observed) else c(0, 1)
      if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)
      circlize::colorRamp2(limits, c("#F1F1F1", "#333333"))
    } else {
      levels <- if (is.factor(v)) levels(v) else sort(unique(as.character(v[!is.na(v)])))
      if (!length(levels)) levels <- "(missing)"
      stats::setNames(grDevices::hcl.colors(length(levels), "Dark 3"), levels)
    }
  })
}

#' Plot a cohort copy-number heatmap
#'
#' Builds an annotated `ComplexHeatmap::Heatmap` from an [ichor_matrix()], one
#' row per sample and one column per genomic bin, split by chromosome.
#'
#' @details
#' Rows keep manifest order unless you cluster or pass `row_order`. Optional
#' clustering uses Euclidean distance and complete linkage on bins observed in
#' every sample, with no imputation; on call codes it assumes ordinal spacing
#' and is exploratory only. Missing, tied and under-covered cells use `na_col`;
#' inspect the matrix `coverage` and `mixed` layers to tell them apart. CN 2 is
#' a reference color, not a neutrality call. See the installed methods
#' contract for limits.
#'
#' `group` splits rows into blocks (for example biofluids or timepoints) while
#' keeping `row_order` within each block; blocks follow factor levels, otherwise
#' first appearance in row order. `row_labels` changes displayed text only: rows
#' stay keyed by sample ID. `text_style` sets font family, size and colour for
#' row labels, block and chromosome titles, annotation names and legends, so a
#' figure can match the text of surrounding panels.
#' @param x An `ichor_matrix`.
#' @param annotation_columns Metadata columns to show as row annotations.
#' @param cluster_rows Explicitly cluster samples? Default `FALSE`.
#' @param show_row_names Display sample identifiers? Default `TRUE`; set
#'   `FALSE` to omit them from the figure.
#' @param colors Custom continuous color function or named discrete call colors.
#' @param annotation_colors Named list of annotation palettes/color functions.
#' @param row_order Optional permutation of sample aliases. Mutually exclusive
#'   with clustering.
#' @param na_col Color for missing, insufficiently covered or tied values.
#' @param group Optional metadata column splitting rows into blocks. No missing
#'   values; cannot be combined with `cluster_rows`.
#' @param row_labels Optional display labels named by every sample ID.
#' @param legend_title Optional title for the value legend.
#' @param legend_direction `"vertical"` (default) or `"horizontal"`, one row.
#' @param text_style Optional list of `grid::gpar()` settings, limited to
#'   `fontfamily`, `fontsize` and `col`.
#' @return An undrawn `ComplexHeatmap::Heatmap` object. Use
#'   `ComplexHeatmap::draw()` to render it on the caller's device.
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'   root <- system.file("extdata", package = "seeNA")
#'   cohort <- read_ichor_cohort(file.path(root, "example-manifest.csv"), "hg38")
#'   m <- ichor_matrix(cohort, value = "call", chromosomes = c("1", "2"))
#'   h <- plot_ichor_heatmap(m, annotation_columns = "condition")
#'   class(h)
#'   # ComplexHeatmap::draw(h) renders it on the current device.
#' }
#' @export
plot_ichor_heatmap <- function(x, annotation_columns = NULL, cluster_rows = FALSE,
                               show_row_names = TRUE, colors = NULL,
                               annotation_colors = NULL, row_order = NULL, na_col = "#D9D9D9",
                               group = NULL, row_labels = NULL, legend_title = NULL,
                               legend_direction = c("vertical", "horizontal"), text_style = NULL) {
  legend_direction <- match.arg(legend_direction)
  validate_ichor_matrix(x)
  .flag(cluster_rows, "cluster_rows")
  .flag(show_row_names, "show_row_names")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    .ichor_abort("Install ComplexHeatmap and circlize with BiocManager to draw heatmaps.", "seena_dependency_error")
  }
  if (!is.null(row_order) && (cluster_rows || anyDuplicated(row_order) ||
                             !setequal(row_order, rownames(x$values)))) .ichor_abort("row_order must be a sample permutation without clustering.")
  rows <- if (is.null(row_order)) seq_len(nrow(x$values)) else match(row_order, rownames(x$values))
  ids <- rownames(x$values)
  labels <- ids
  if (!is.null(row_labels)) {
    if (!is.character(row_labels) || is.null(names(row_labels)) || anyNA(row_labels) ||
        !all(ids %in% names(row_labels))) .ichor_abort("row_labels must be a character vector named by every sample ID.")
    labels <- unname(row_labels[ids])
  }
  split <- NULL
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || !group %in% setdiff(names(x$samples), "sample_id")) {
      .ichor_abort("group must name one metadata column of the matrix samples.")
    }
    if (cluster_rows) .ichor_abort("group cannot be combined with cluster_rows; supply row_order instead.")
    g <- x$samples[[group]][match(ids, x$samples$sample_id)]
    if (anyNA(g)) .ichor_abort("group column has missing values.")
    split <- factor(as.character(g), levels = if (is.factor(g)) levels(droplevels(g)) else unique(as.character(g[rows])))
  }
  if (!is.null(legend_title) && (!is.character(legend_title) || length(legend_title) != 1L || is.na(legend_title))) {
    .ichor_abort("legend_title must be one string.")
  }
  if (is.null(text_style)) text_style <- list()
  if (!is.list(text_style) || (length(text_style) && (is.null(names(text_style)) ||
      !all(names(text_style) %in% c("fontfamily", "fontsize", "col"))))) {
    .ichor_abort("text_style may only set fontfamily, fontsize and col.")
  }
  gp <- function(...) do.call(grid::gpar, utils::modifyList(list(...), text_style))
  clustering <- FALSE
  if (cluster_rows) {
    if (nrow(x$values) < 2 || nrow(x$values) > 2000) .ichor_abort("Clustering requires 2-2000 samples; supply a row order otherwise.")
    complete <- colSums(is.na(x$values)) == 0
    if (sum(complete) < 2) .ichor_abort("Clustering needs at least two bins observed in every sample; no imputation is performed.")
    distances <- stats::dist(x$values[, complete, drop = FALSE])
    if (any(!is.finite(distances))) .ichor_abort("Non-finite clustering distances.")
    clustering <- stats::hclust(distances, method = "complete")
  }
  if (is.null(colors)) {
    observed <- x$values[is.finite(x$values)]
    if (x$value == "call") {
      colors <- c("-2" = "#A6E4C0", "-1" = "#1E9E5E", "0" = "#F7F7F7", "1" = "#CB2B2B")
    } else if (x$value %in% c("copy_number", "corrected_copy_number")) {
      colors <- circlize::colorRamp2(c(0, 2, max(4, observed)), c("#1E9E5E", "#F7F7F7", "#CB2B2B"))
    } else {
      extent <- max(1, abs(observed))
      colors <- circlize::colorRamp2(c(-extent, 0, extent), c("#1E9E5E", "#F7F7F7", "#CB2B2B"))
    }
  }
  if (x$value == "call") .check_colors(colors, c("-2", "-1", "0", "1"))
  left <- NULL
  if (length(annotation_columns)) {
    if (!all(annotation_columns %in% names(x$samples))) .ichor_abort("Unknown annotation columns.")
    ann <- x$samples[match(rownames(x$values), x$samples$sample_id), annotation_columns, drop = FALSE]
    rownames(ann) <- rownames(x$values)
    palettes <- .annotation_palette(ann)
    if (!is.null(annotation_colors)) {
      if (!is.list(annotation_colors) || is.null(names(annotation_colors)) ||
          !all(names(annotation_colors) %in% names(ann))) .ichor_abort("annotation_colors must be named by annotation column.")
      palettes[names(annotation_colors)] <- annotation_colors
    }
    left <- ComplexHeatmap::rowAnnotation(df = ann, col = palettes, na_col = na_col, annotation_name_gp = gp(),
      annotation_legend_param = lapply(ann, function(v) list(title_gp = gp(fontface = "bold"), labels_gp = gp())))
  }
  legend <- if (x$value == "call") {
    list(title = paste(x$call_column, x$genome_build, sep = "\n"),
         at = c("-2", "-1", "0", "1"), labels = c("Deep loss", "Loss", "Neutral", "Gain"))
  } else list(title = paste(x$value, x$genome_build, sep = "\n"))
  if (!is.null(legend_title)) legend$title <- legend_title
  if (legend_direction == "horizontal") legend <- c(legend, list(direction = "horizontal", nrow = 1))
  # outlined keys: a white or near-white call state is otherwise invisible
  if (x$value == "call") legend$border <- "grey60"
  legend$title_gp <- gp(fontface = "bold")
  legend$labels_gp <- gp()
  raster <- ncol(x$values) > 1000 && requireNamespace("ragg", quietly = TRUE)
  suppressMessages(ComplexHeatmap::Heatmap(
    x$values, name = x$value, col = colors, na_col = na_col,
    cluster_rows = clustering, row_order = if (cluster_rows) NULL else rows,
    row_split = split, cluster_row_slices = FALSE, row_title_gp = gp(fontface = "bold"),
    row_labels = labels, row_names_gp = gp(),
    cluster_columns = FALSE, cluster_column_slices = FALSE,
    column_split = factor(x$bins$chr, levels = .chr_levels),
    show_column_names = FALSE, show_row_names = show_row_names,
    use_raster = raster, raster_device = if (raster) "agg_png" else "png",
    raster_resize_mat = FALSE, raster_by_magick = FALSE,
    left_annotation = left, heatmap_legend_param = legend,
    column_title = unique(x$bins$chr), column_title_side = "bottom",
    column_title_gp = gp(fontsize = 9), border = FALSE
  ))
}
