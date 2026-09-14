test_that("profile, comparison and region plots render as expected", {
  skip_if_not_installed("vdiffr")
  a <- example_sample("a")
  b <- example_sample("b")
  vdiffr::expect_doppelganger("profile", plot_ichor_profile(a))
  vdiffr::expect_doppelganger("profile-ploidy-adjusted", plot_ichor_profile(a, ploidy_adjust = TRUE))
  vdiffr::expect_doppelganger("compare", plot_ichor_compare(list(a, b)))
  vdiffr::expect_doppelganger("region", plot_ichor_region(list(a, b), "chr1:1-4000000"))
})
