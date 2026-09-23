window_file <- function(chr = "19", terminal_start = 58000001, terminal_end = 59000000) {
  d <- data.table::fread(fixture_path("example-a.cna.seg"), data.table = FALSE)[1:3, ]
  d$chr <- chr
  d$start <- c(1, 1000001, terminal_start)
  d$end <- c(1000000, 2000000, terminal_end)
  f <- tempfile(fileext = ".cna.seg")
  data.table::fwrite(d, f, sep = "\t")
  f
}

test_that("window policy recognizes regular terminal padding and records it", {
  f <- window_file()
  expect_message(s <- read_ichor_sample(f, genome_build = "hg38"), "Trimmed 1")
  expect_identical(s$coordinate_policy, "window")
  expect_equal(s$bins$end[3], 58617616)
  expect_equal(s$coordinate_changes$original_end, 59000000)
  expect_error(read_ichor_sample(f, genome_build = "hg38", bounds = "error"), "out-of-build")
  expect_message(m <- read_ichor_cohort(data.frame(sample_id = "alias", cna_seg = f), "hg38"), "across 1 sample")
  expect_equal(ichor_matrix(m, value = "logR", chromosomes = "19")$coverage[1, 59], 1)
})

test_that("window policy refuses off-grid, excessive or unsupported padding", {
  expect_error(read_ichor_sample(window_file(terminal_end = 59000001), genome_build = "hg38"), class = "seena_bounds_error")
  expect_error(read_ichor_sample(window_file(terminal_start = 58000002), genome_build = "hg38"), class = "seena_bounds_error")
  expect_error(read_ichor_sample(window_file(terminal_end = 60000000), genome_build = "hg38"), class = "seena_bounds_error")
  # Both builds can share the same terminal window: acceptance is not build proof.
  f <- window_file(chr = "3", terminal_start = 198000001, terminal_end = 199000000)
  expect_message(a <- read_ichor_sample(f, genome_build = "hg19"), "caller-specified")
  expect_message(b <- read_ichor_sample(f, genome_build = "hg38"), "caller-specified")
  expect_false(identical(a$bins$end, b$bins$end))
})

test_that("cohort aggregates only its own coordinate notices", {
  f <- window_file()
  tab <- data.frame(sample_id = c("a", "b"), cna_seg = c(f, f))
  messages <- character(); warnings <- character()
  c <- withCallingHandlers(read_ichor_cohort(tab, "hg38"),
    message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_length(messages, 1)
  expect_match(messages, "Trimmed 2.*across 2")
  expect_length(c$samples, 2)
  withCallingHandlers(read_ichor_cohort(tab, "hg38", bounds = "trim"),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_length(warnings, 1)
  expect_match(warnings, "Trimmed 2")
})

test_that("segment padding must match an observed terminal bin boundary", {
  f <- window_file()
  d <- data.table::fread(fixture_path("example-a.seg"), data.table = FALSE)[1, ]
  d$chr <- "19"; d$start <- 1; d$end <- 59000000
  seg <- tempfile(fileext = ".seg")
  data.table::fwrite(d, seg, sep = "\t")
  expect_message(s <- read_ichor_sample(f, seg, genome_build = "hg38"), "Trimmed 2")
  expect_equal(s$segments$end, s$bins$end[3])
  expect_setequal(s$coordinate_changes$role, c("bins", "segments"))
  d$end <- 59000001
  data.table::fwrite(d, seg, sep = "\t")
  expect_error(read_ichor_sample(f, seg, genome_build = "hg38"), class = "seena_bounds_error")
})

test_that("cohort padding aggregation does not suppress unrelated warnings", {
  fingerprint <- .fingerprint
  local_mocked_bindings(.fingerprint = function(paths) {
    warning("unrelated import warning", call. = FALSE)
    fingerprint(paths)
  })
  tab <- data.frame(sample_id = "alias", cna_seg = window_file())
  seen <- character()
  expect_message(withCallingHandlers(read_ichor_cohort(tab, "hg38"),
    warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") }), "Trimmed 1")
  expect_identical(seen, rep("unrelated import warning", 2))
})

test_that("parameter accessors and adjusted logR are explicit and immutable", {
  a <- example_sample("a")
  original <- a
  expect_equal(ichor_tf(a), 0.125)
  expect_equal(ichor_ploidy(a), 2.1)
  shift <- log2((0.125 * 2.1 + (1 - 0.125) * 2) / 2)
  expect_equal(ichor_adjusted_logr(a), a$bins$logR + shift)
  expect_equal(ichor_adjusted_logr(a, "segments"), a$segments$median + shift)
  expect_identical(a, original)
  a$bins$logR[1] <- NA_real_
  expect_true(is.na(ichor_adjusted_logr(a)[1]))
  a$params$tumor_fraction <- 0
  expect_equal(ichor_adjusted_logr(a), a$bins$logR)
  a$params <- NULL
  expect_identical(ichor_tf(a), NA_real_)
  expect_identical(ichor_ploidy(a), NA_real_)
  expect_error(ichor_adjusted_logr(a), class = "seena_parameter_error")
  expect_error(plot_ichor_profile(a, ploidy_adjust = TRUE), class = "seena_parameter_error")
})

test_that("all plotting paths apply the same optional shift to bins and segments", {
  a <- example_sample("a"); b <- example_sample("b")
  p <- plot_ichor_profile(a, ploidy_adjust = TRUE)
  expect_equal(p$data$logR, ichor_adjusted_logr(a))
  expect_equal(p$layers[[4]]$data$median, ichor_adjusted_logr(a, "segments"))
  compare <- plot_ichor_compare(list(a, b), ploidy_adjust = TRUE)
  expect_equal(compare$data$logR[compare$data$sample == a$sample_id], ichor_adjusted_logr(a))
  seg <- compare$layers[[3]]$data
  expect_equal(seg$median[seg$sample == b$sample_id], ichor_adjusted_logr(b, "segments"))
  region <- plot_ichor_region(a, "chr1:1-4000000", ploidy_adjust = TRUE)
  expect_equal(region$data$logR, ichor_adjusted_logr(a)[1:4])
  expect_true(attr(region, "ichor_transform")$ploidy_adjust)
  expect_false(attr(plot_ichor_profile(a), "ichor_transform")$ploidy_adjust)
})

test_that("state colors and genome layout share a public contract", {
  expect_identical(levels(ichor_call_state(c("AMP", "HLAMP3", NA))), names(ichor_state_colors()))
  d <- ichor_genome_layout("hg38")
  expect_identical(d$chr, c(as.character(1:22), "X"))
  expect_equal(d$offset, c(0, head(cumsum(d$length), -1)))
  expect_equal(d$mid, d$offset + d$length / 2)
  expect_identical(ichor_genome_layout("hg19", c("chrY", "1", "X"))$chr, c("1", "X", "Y"))
  expect_error(ichor_genome_layout("hg38", "MT"), "Unsupported")
})

test_that("neutral evidence never guesses absent or conflicting baselines", {
  a <- example_sample("a")
  n <- ichor_neutral_cn(a)
  expect_equal(n$neutral_cn, c(2, 2, NA))
  expect_identical(n$status, c("supported", "supported", "no_neutral_bins"))
  # Sex metadata do not change the observed CN reference.
  a$params$gender <- "male"
  expect_equal(ichor_neutral_cn(a)$neutral_cn, n$neutral_cn)
  xneutral <- which(a$bins$chr == "X" & a$bins$corrected_call == "NEUT")
  a$bins$corrected_copy_number[xneutral] <- 1
  expect_equal(ichor_neutral_cn(a)$neutral_cn, c(2, 1, NA))
  a$bins$corrected_copy_number[xneutral[1]] <- 2
  n <- ichor_neutral_cn(a)
  expect_true(is.na(n$neutral_cn[2]))
  expect_identical(n$status[2], "ambiguous")
  expect_equal(n$candidates[[2]], c(1, 2))
  a$bins$corrected_call[a$bins$chr == "X"] <- "GAIN"
  expect_identical(ichor_neutral_cn(a)$status[2], "no_neutral_bins")
  a$bins$corrected_copy_number[a$bins$chr != "X"] <- NA_real_
  expect_identical(ichor_neutral_cn(a)$status[1], "missing_cn")
  expect_equal(ichor_neutral_cn(a, call_column = "event")$neutral_cn, c(2, 2, NA))
})
