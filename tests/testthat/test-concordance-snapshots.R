test_that("paired concordance profiles and unknown evidence render consistently", {
  skip_if_not_installed("vdiffr")
  a <- example_sample("a"); b <- example_sample("b")
  vdiffr::expect_doppelganger("concordance-region",
    plot_ichor_concordance(a, b, region = "chr1:1-5000000", point_size = 1.5,
                           ploidy_adjust = TRUE, show_sample_id = FALSE))
  vdiffr::expect_doppelganger("concordance-styled",
    plot_ichor_concordance(a, b, region = "chr1:1-5000000", point_size = 1.5,
      heights = c(2, 1.1), segment_linewidth = 0.32, point_stroke = 0,
      show_sample_id = FALSE) + ggplot2::theme(legend.position = "bottom"))
  # Modified format fixtures exercise opposite-direction bars and X evidence,
  # not a patient example or a new ichorCNA fit.
  a$bins$corrected_call[9:12] <- c("HETD", "GAIN", "NEUT", "GAIN")
  b$bins$corrected_call[9:12] <- "GAIN"
  a$bins$logR[9:12] <- c(-0.5, 0.5, 0, 0)
  b$bins$logR[9:12] <- c(0.5, 0.5, 0.5, 0.5)
  vdiffr::expect_doppelganger("concordance-sex-reference",
    plot_ichor_concordance(a, b, region = "chrX:1-5000000", point_size = 1.5,
                           show_sample_id = FALSE))
})
