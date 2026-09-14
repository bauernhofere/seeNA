test_that("numeric-looking file-manifest aliases retain leading zeroes", {
  f <- tempfile(fileext = ".csv")
  tab <- data.frame(sample_id = c("001", "002"),
                    cna_seg = c(fixture_path("example-a.cna.seg"), fixture_path("example-b.cna.seg")))
  data.table::fwrite(tab, f)
  c <- read_ichor_cohort(f, "hg38")
  expect_identical(c$metadata$sample_id, c("001", "002"))
})

test_that("matrix validators reject scientific metadata drift", {
  c <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  original <- ichor_matrix(c, value = "call", chromosomes = "1")
  m <- original; m$values[1, 1] <- 0.5
  expect_error(validate_ichor_matrix(m), "Fractional")
  m <- original; m$aggregation <- "mean"
  expect_error(validate_ichor_matrix(m), "aggregation")
  m <- original; m$coverage[1, 1] <- 0.5
  expect_error(validate_ichor_matrix(m), "coverage")
  m <- original; m$bins$end[1] <- m$bins$end[1] - 1
  expect_error(validate_ichor_matrix(m), "target grid")
  m <- original; m$call_column <- "auto"
  expect_error(validate_ichor_matrix(m), "call column")
})

test_that("normalization and hg19 coordinates are not accidentally hg38-specific", {
  expect_equal(.normalize_chr(c("chr1", "23", "chr24", NA)), c("1", "X", "Y", NA))
  expect_equal(parse_ichor_region("chr1", "hg19")$end, 249250621)
  expect_equal(parse_ichor_region("chrY", "hg38")$end, 57227415)
  expect_error(parse_ichor_region("chr1:248956423-249000000", "hg38"), "outside")
  expect_equal(parse_ichor_region("chr1:248956423-249000000", "hg19")$start, 248956423)
})

test_that("unknown state rendering cannot borrow the neutral color", {
  a <- example_sample("a")
  a$bins$corrected_call[1] <- NA_character_
  p <- plot_ichor_profile(a)
  built <- ggplot2::ggplot_build(p)
  points <- built$data[[3]]
  expect_equal(points$colour[1], "grey65")
  expect_false(points$colour[1] == ichor_state_colors()[["Neutral"]])
  expect_equal(p$scales$get_scales("x")$limits[2],
               sum(.chromosome_sizes("hg38")$length[1:23]))
})

test_that("a sample with Y segments but no Y bins can be drawn", {
  a <- example_sample("a")
  last <- a$segments[nrow(a$segments), ]
  last$chr <- "Y"
  a$segments <- rbind(a$segments, last)
  expect_invisible(validate_ichor_sample(a))
  expect_s3_class(ggplot2::ggplot_build(plot_ichor_profile(a)), "ggplot_built")
})

test_that("ggplot and rasterized heatmaps render headlessly", {
  skip_if_not_installed("ragg")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  c <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  m <- ichor_matrix(c, value = "call")
  f <- tempfile(fileext = ".png")
  ragg::agg_png(f, width = 1200, height = 600)
  tryCatch(ComplexHeatmap::draw(plot_ichor_heatmap(m, annotation_columns = "condition")),
           finally = grDevices::dev.off())
  expect_gt(file.info(f)$size, 1000)
  g <- tempfile(fileext = ".png")
  ggplot2::ggsave(g, plot_ichor_profile(c$samples[[1]]), width = 8, height = 3,
                 device = ragg::agg_png)
  expect_gt(file.info(g)$size, 1000)
})
