# Cohort heatmap -----------------------------------------------------------

#' Plot a cohort copy-number heatmap
#'
#' Uses ComplexHeatmap when installed. Samples are rows and ordered genomic bins
#' are columns split by chromosome.
#'
#' @param x An `ichor_matrix`.
#' @param annotation_columns Optional columns from the cohort metadata to show
#'   as row annotations.
#' @param cluster_rows Whether to cluster samples.
#' @param show_row_names Whether to display sample identifiers.
#' @param colors Optional color function. Defaults depend on the matrix value.
#' @return A `ComplexHeatmap::Heatmap` object.
#' @export
plot_ichor_heatmap <- function(x, annotation_columns = NULL, cluster_rows = TRUE,
                               show_row_names = nrow(x$values) <= 100, colors = NULL) {
  if (!inherits(x, "ichor_matrix")) .ichor_abort("x must be an ichor_matrix.")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    .ichor_abort("Install ComplexHeatmap and circlize to draw cohort heatmaps.",
                 "ichorviz_dependency_error")
  }
  if (is.null(colors)) {
    colors <- if (x$value == "call") {
      circlize::colorRamp2(c(-2, -1, 0, 1), c("#A6E4C0", "#1E9E5E", "#F7F7F7", "#CB2B2B"))
    } else if (x$value %in% c("copy_number", "corrected_copy_number")) {
      circlize::colorRamp2(c(0, 2, 4), c("#1E9E5E", "#F7F7F7", "#CB2B2B"))
    } else {
      circlize::colorRamp2(c(-1, 0, 1), c("#1E9E5E", "#F7F7F7", "#CB2B2B"))
    }
  }
  left <- NULL
  if (!is.null(annotation_columns)) {
    missing <- setdiff(annotation_columns, names(x$samples))
    if (length(missing)) .ichor_abort(sprintf("Unknown annotation column%s: %s",
                                               if (length(missing) > 1) "s" else "",
                                               paste(missing, collapse = ", ")))
    ann <- x$samples[annotation_columns]
    rownames(ann) <- x$samples$sample_id
    left <- ComplexHeatmap::rowAnnotation(df = ann)
  }
  suppressMessages(ComplexHeatmap::Heatmap(
    x$values,
    name = x$value,
    col = colors,
    cluster_rows = cluster_rows,
    cluster_columns = FALSE,
    column_split = factor(x$bins$chr, levels = .chr_levels),
    show_column_names = FALSE,
    show_row_names = show_row_names,
    use_raster = ncol(x$values) > 1000,
    raster_device = if (requireNamespace("ragg", quietly = TRUE)) "agg_png" else "png",
    left_annotation = left,
    column_title = NULL,
    border = FALSE
  ))
}
