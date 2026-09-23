example_call_matrix <- function() {
  cohort <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  ichor_matrix(cohort, value = "call")
}

drawn_row_order <- function(h) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  ComplexHeatmap::row_order(ComplexHeatmap::draw(h))
}

test_that("group splits rows into blocks in first-appearance or factor order", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  m <- example_call_matrix()
  ids <- m$samples$sample_id

  h <- plot_ichor_heatmap(m, group = "condition", row_order = rev(ids))
  expect_identical(drawn_row_order(h), list(B = 2L, A = 1L))

  m$samples$condition <- factor(m$samples$condition, levels = c("A", "B"))
  h <- plot_ichor_heatmap(m, group = "condition", row_order = rev(ids))
  expect_identical(drawn_row_order(h), list(A = 1L, B = 2L))
})

test_that("row labels change text only and text style reaches every text element", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  m <- example_call_matrix()
  style <- list(fontfamily = "sans", fontsize = 7, col = "grey20")
  h <- plot_ichor_heatmap(m, row_labels = c("example-b" = "same", "example-a" = "same"),
                          legend_title = "Call", text_style = style, annotation_columns = "condition")
  expect_identical(h@row_names_param$labels, c("same", "same"))
  expect_identical(nrow(h@matrix), 2L)
  expect_identical(h@matrix_legend_param$title, "Call")
  for (g in list(h@row_names_param$gp, h@matrix_legend_param$title_gp, h@matrix_legend_param$labels_gp)) {
    expect_identical(g$fontfamily, "sans")
    expect_equal(g$fontsize, 7)
    expect_identical(g$col, "grey20")
  }
  expect_equal(unname(h@matrix_legend_param$title_gp$font), 2L)
  expect_identical(h@matrix_legend_param$border, "grey60")
  expect_null(h@matrix_legend_param$direction)
  horizontal <- plot_ichor_heatmap(m, legend_direction = "horizontal")
  expect_identical(horizontal@matrix_legend_param$direction, "horizontal")
  expect_equal(horizontal@matrix_legend_param$nrow, 1)

  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  tryCatch(ComplexHeatmap::draw(plot_ichor_heatmap(m, group = "condition", text_style = style)),
           finally = grDevices::dev.off())
  expect_gt(file.info(f)$size, 1000)
})

test_that("grouping, labels and style arguments are validated", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  m <- example_call_matrix()
  expect_error(plot_ichor_heatmap(m, group = "missing"), "group")
  expect_error(plot_ichor_heatmap(m, group = "sample_id"), "group")
  expect_error(plot_ichor_heatmap(m, group = "condition", cluster_rows = TRUE), "cluster_rows")
  bad <- m
  bad$samples$condition[1] <- NA
  expect_error(plot_ichor_heatmap(bad, group = "condition"), "missing values")
  expect_error(plot_ichor_heatmap(m, row_labels = c("example-a" = "x")), "row_labels")
  expect_error(plot_ichor_heatmap(m, legend_title = c("a", "b")), "legend_title")
  expect_error(plot_ichor_heatmap(m, legend_direction = "diagonal"))
  expect_error(plot_ichor_heatmap(m, text_style = list(fontface = "bold")), "text_style")
  expect_error(plot_ichor_heatmap(m, text_style = list(7)), "text_style")
})
