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
#' Default row order is manifest order. Optional clustering uses Euclidean
#' distances and complete linkage on columns observed in every sample (at least
#' two required), without imputation. Calls use discrete colors; continuous
#' default scales include the observed range. CN 2 is a visual reference, not a
#' neutrality call; sex-chromosome baseline and fitted ploidy may differ.
#' @param x An `ichor_matrix`.
#' @param annotation_columns Metadata columns to show as row annotations.
#' @param cluster_rows Explicitly cluster samples? Default FALSE.
#' @param show_row_names Display sample aliases? Default FALSE.
#' @param colors Custom continuous color function or named discrete call colors.
#' @param annotation_colors Named list of annotation palettes/color functions.
#' @param row_order Optional permutation of sample aliases. Mutually exclusive
#'   with clustering.
#' @param na_col Color for missing, insufficiently covered or tied values.
#' @return An undrawn `ComplexHeatmap::Heatmap` object. Use
#'   `ComplexHeatmap::draw()` to render it on the caller's device.
#' @export
plot_ichor_heatmap <- function(x, annotation_columns = NULL, cluster_rows = FALSE,
                               show_row_names = FALSE, colors = NULL,
                               annotation_colors = NULL, row_order = NULL, na_col = "#D9D9D9") {
  validate_ichor_matrix(x)
  .flag(cluster_rows, "cluster_rows")
  .flag(show_row_names, "show_row_names")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    .ichor_abort("Install ComplexHeatmap and circlize with BiocManager to draw heatmaps.", "ichorviz_dependency_error")
  }
  if (!is.null(row_order) && (cluster_rows || anyDuplicated(row_order) ||
                             !setequal(row_order, rownames(x$values)))) .ichor_abort("row_order must be a sample permutation without clustering.")
  rows <- if (is.null(row_order)) seq_len(nrow(x$values)) else match(row_order, rownames(x$values))
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
    left <- ComplexHeatmap::rowAnnotation(df = ann, col = palettes, na_col = na_col)
  }
  legend <- if (x$value == "call") {
    list(title = paste(x$call_column, x$genome_build, sep = "\n"),
         at = c("-2", "-1", "0", "1"), labels = c("Deep loss", "Loss", "Neutral", "Gain"))
  } else list(title = paste(x$value, x$genome_build, sep = "\n"))
  raster <- ncol(x$values) > 1000 && requireNamespace("ragg", quietly = TRUE)
  suppressMessages(ComplexHeatmap::Heatmap(
    x$values, name = x$value, col = colors, na_col = na_col,
    cluster_rows = clustering, row_order = if (cluster_rows) NULL else rows,
    cluster_columns = FALSE, cluster_column_slices = FALSE,
    column_split = factor(x$bins$chr, levels = .chr_levels),
    show_column_names = FALSE, show_row_names = show_row_names,
    use_raster = raster, raster_device = if (raster) "agg_png" else "png",
    raster_resize_mat = FALSE, raster_by_magick = FALSE,
    left_annotation = left, heatmap_legend_param = legend,
    column_title = unique(x$bins$chr), column_title_side = "bottom",
    column_title_gp = grid::gpar(fontsize = 9), border = FALSE
  ))
}
