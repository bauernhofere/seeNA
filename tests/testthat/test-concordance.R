test_that("pair classifications follow calls, not logR sign or CN magnitude", {
  a <- example_sample("a"); b <- example_sample("b")
  a$bins$corrected_call[1:6] <- c("NEUT", "GAIN", "NEUT", "HOMD", "GAIN", NA)
  b$bins$corrected_call[1:6] <- c("NEUT", "NEUT", "HETD", "HETD", "HOMD", "GAIN")
  a$bins$logR[5] <- -0.8; b$bins$logR[5] <- 0.2
  d <- ichor_pair_concordance(a, b, chromosomes = c("1", "2"))
  observed <- d[d$coverage_a > 0 | d$coverage_b > 0, ]
  expect_identical(observed$concordance[1:6], c("both_neutral", "a_only", "b_only", "concordant", "discordant", "unknown"))
  expect_equal(observed$representative_logR[5], -0.8)
  expect_identical(observed$representative_source[5], "a")
  expect_true(any(d$coverage_a == 0 & d$concordance == "unknown"))
  expect_identical(attr(d, "settings")$call_column, "corrected_call")
  expect_identical(attr(d, "provenance"), ichor_provenance(structure(list(samples = list(a, b)), class = "ichor_cohort")))
})

test_that("pair alignment handles shifted bins with coverage and mixed-call diagnostics", {
  a <- example_sample("a"); b <- example_sample("b")
  b$segments <- NULL
  i <- which(b$bins$chr == "1")
  b$bins$start[i] <- b$bins$start[i] + 5e5
  b$bins$end[i] <- b$bins$end[i] + 5e5
  b$bins$corrected_call[i] <- c("NEUT", "GAIN", "GAIN", "GAIN")
  b$bins$logR[i] <- c(0, 1, 1, 1)
  d <- ichor_pair_concordance(a, b, chromosomes = "1")
  expect_equal(d$coverage_b[1:5], c(0.5, 1, 1, 1, 0.5))
  expect_identical(d$concordance[1], "unknown")
  expect_true(d$mixed_b[2])
  expect_identical(d$concordance[2], "unknown")
  expect_match(d$reason[2], "mixed_calls_b")
  expect_equal(d$raw_logR_b[2:4], c(0.5, 1, 1))
  expect_identical(d$concordance[3:4], c("concordant", "concordant"))
  relaxed <- ichor_pair_concordance(a, b, chromosomes = "1", min_coverage = 0.5)
  expect_identical(relaxed$concordance[1], "both_neutral")
})

test_that("disjoint partial support is not a co-observed comparison", {
  a <- example_sample("a"); b <- example_sample("b")
  a$bins <- a$bins[1, , drop = FALSE]; b$bins <- b$bins[1, , drop = FALSE]
  a$segments <- b$segments <- NULL
  a$bins$end <- 500000; b$bins$start <- 500001
  a$bins$corrected_call <- "GAIN"
  d <- ichor_pair_concordance(a, b, chromosomes = "1", min_coverage = 0.5)
  expect_equal(c(d$coverage_a[1], d$coverage_b[1]), c(0.5, 0.5))
  expect_equal(d$joint_coverage[1], 0)
  expect_identical(d$call_a[1], "Gain")
  expect_identical(d$concordance[1], "unknown")
  expect_match(d$reason[1], "insufficient_joint_call_support")
})

test_that("unknown mixed calls cannot be rescued by a winning mode", {
  a <- example_sample("a"); b <- example_sample("b")
  a$bins$end[1] <- 750000
  a$bins$start[2] <- 750001
  a$bins$end[2] <- 1000000
  a$bins$corrected_call[1:2] <- c("GAIN", "NEUT")
  d <- ichor_pair_concordance(a, b, chromosomes = "1")
  expect_identical(d$call_a[1], "Gain")
  expect_equal(d$coverage_a[1], 1)
  expect_true(d$mixed_a[1])
  expect_identical(d$concordance[1], "unknown")
})

test_that("missing logR does not erase a known call classification", {
  a <- example_sample("a"); b <- example_sample("b")
  a$bins$logR[3] <- NA_real_
  d <- ichor_pair_concordance(a, b, chromosomes = "1")
  expect_identical(d$concordance[3], "a_only")
  expect_equal(d$coverage_a[3], 1)
  expect_equal(d$logR_coverage_a[3], 0)
  expect_true(is.na(d$representative_logR[3]))
  p <- plot_ichor_concordance(a, b, region = "chr1:1-4000000")
  expect_true(any(p$layers[[1]]$data$xleft < 2.1 & p$layers[[1]]$data$xright > 2.9))
})

test_that("both heights and representative rule expose discordant ties", {
  a <- example_sample("a"); b <- example_sample("b")
  a$bins$corrected_call[1:2] <- "GAIN"
  b$bins$corrected_call[1:2] <- c("HETD", "GAIN")
  a$bins$logR[1:2] <- c(0.5, 0.4); b$bins$logR[1:2] <- c(-0.5, 0.8)
  d <- ichor_pair_concordance(a, b, chromosomes = "1")
  swapped <- ichor_pair_concordance(b, a, chromosomes = "1")
  expect_equal(d$representative_logR[1:2], c(0.5, 0.6))
  expect_equal(swapped$representative_logR[1:2], c(-0.5, 0.6))
  expect_equal(d$representative_source[1:2], c("a", "mean"))
  p <- plot_ichor_concordance(a, b, region = "chr1:1-2000000")
  bars <- p$layers[[4]]$data
  expect_equal(sort(bars$value), sort(c(0.5, -0.5, 0.4, 0.8)))
  expect_identical(attr(p, "ichor_view")$height, "both")
  rep <- plot_ichor_concordance(a, b, region = "chr1:1-2000000", height = "representative")
  expect_equal(rep$layers[[4]]$data$value, c(0.5, 0.6))
})

test_that("sex reference absence is flagged rather than disproving source calls", {
  a <- example_sample("a"); b <- example_sample("b")
  b$bins$corrected_call[b$bins$chr == "X"] <- "GAIN"
  d <- ichor_pair_concordance(a, b, chromosomes = "X")
  expect_true(all(d$baseline_flag))
  expect_identical(d$baseline_status_b[1], "no_neutral_bins")
  expect_identical(d$concordance[1:4], c("b_only", "concordant", "concordant", "b_only"))
  strict <- ichor_pair_concordance(a, b, chromosomes = "X", sex_chromosomes = "require_neutral")
  expect_true(all(strict$concordance == "unknown"))
  expect_identical(strict$call_b[1], "Gain")
  expect_match(strict$reason[1], "sex_reference_unresolved")
  p <- plot_ichor_concordance(a, b, region = "chrX:1-4000000")
  expect_equal(nrow(p$layers[[5]]$data), 4)
  b$bins$corrected_copy_number <- NULL
  expect_identical(ichor_pair_concordance(a, b, chromosomes = "X")$baseline_status_b[1], "unavailable_columns")
})

test_that("pair shifts, identities, parameters and limits are explicit", {
  a <- example_sample("a"); b <- example_sample("b")
  original <- a
  d <- ichor_pair_concordance(a, b, chromosomes = "1", ploidy_adjust = TRUE)
  expect_equal(d$logR_a[1:4], ichor_adjusted_logr(a)[1:4])
  expect_identical(a, original)
  expect_error(ichor_pair_concordance(a, a), "unique")
  b$genome_build <- "hg19"
  expect_error(ichor_pair_concordance(a, b), "same genome")
  b <- example_sample("b"); b$params <- NULL
  expect_error(ichor_pair_concordance(a, b, ploidy_adjust = TRUE), "requires finite")
  expect_error(ichor_pair_concordance(a, b, max_cells = 1), "max_cells")
  expect_error(plot_ichor_concordance(a, b, ylim = c(0, 1)), "straddle")
  expect_warning(p <- plot_ichor_concordance(a, b, region = "chr1:1-4000000", ylim = c(-0.1, 0.1)), "clipped")
  expect_gt(attr(p, "ichor_view")$clipped_heights, 0)
  expect_gt(max(attr(p, "ichor_concordance")$logR_a), 0.1)
})

test_that("concordance plot uses original upper midpoints and separate legends", {
  a <- example_sample("a"); b <- example_sample("b")
  p <- plot_ichor_concordance(a, b, region = "chr1:750000-1500000", show_sample_id = FALSE)
  expect_equal(p$layers[[3]]$data$x, rep(c(0.5000005, 1.5000005), 2))
  expect_false(grepl("example-", p$labels$subtitle))
  expect_match(p$labels$subtitle, "TF")
  expect_true(all(c("colour", "fill") %in% unlist(lapply(p$scales$scales, `[[`, "aesthetics"))))
  built <- ggplot2::ggplot_build(p)
  expect_identical(as.character(built$layout$layout$panel), c("Profiles", "Directional call agreement"))
})

test_that("zero altered heights have a marker rather than disappearing as neutral", {
  a <- example_sample("a"); b <- example_sample("b")
  a$bins$logR[3] <- 0
  p <- plot_ichor_concordance(a, b, region = "chr1:1-4000000")
  expect_equal(p$layers[[6]]$data$value, 0)
  expect_identical(p$layers[[6]]$data$category, "a_only")
})

test_that("concordance track renders headlessly", {
  a <- example_sample("a"); b <- example_sample("b")
  p <- plot_ichor_concordance(a, b, region = "chr1:1-5000000")
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f, width = 10, height = 6)
  tryCatch(print(p), finally = grDevices::dev.off())
  expect_gt(file.info(f)$size, 0)
})
