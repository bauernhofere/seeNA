test_that("missing is never neutral or false and unknown tokens fail", {
  expect_equal(.call_score(c(NA, "NEUT", "HETD", "GAIN", "HLAMP12", "HOMD")), c(NA, 0, -1, 1, 1, -2))
  expect_equal(.as_flag(c(NA, "", "false", "TRUE", "0", "1")), c(NA, NA, FALSE, TRUE, FALSE, TRUE))
  expect_error(ichor_call_state("TYPO"), class = "ichorviz_schema_error")
  expect_error(.as_flag("uncertain"), class = "ichorviz_schema_error")
  expect_error(.as_number("bad"), class = "ichorviz_schema_error")
  expect_error(.as_number("Inf"), class = "ichorviz_schema_error")
  expect_equal(.as_number(c("NA", "1.5", "")), c(NA, 1.5, NA))
})

test_that("source identity cannot be overridden by an alias", {
  expect_error(read_ichor_sample(fixture_path("example-a.cna.seg"),
    fixture_path("example-b.seg"), fixture_path("example-b.params.txt"),
    "hg38", sample_id = "alias"), class = "ichorviz_identity_error")
  a <- read_ichor_sample(fixture_path("example-a.cna.seg"), fixture_path("example-a.seg"),
    fixture_path("example-a.params.txt"), "hg38", sample_id = "alias")
  expect_identical(a$sample_id, "alias")
  expect_identical(a$params$sample_id, "alias")
  expect_null(a$segments$sample_id)
  expect_null(attr(a$bins, "source_id"))
  expect_false("path" %in% names(ichor_provenance(a)))
  # NULL alone means "use the source ID"; an NA alias is an explicit error.
  expect_error(read_ichor_sample(fixture_path("example-a.cna.seg"), genome_build = "hg38",
    sample_id = NA_character_), "non-empty character scalar")
})

test_that("all interval layers and measurements are validated", {
  for (layer in c("bins", "segments")) {
    a <- example_sample("a")
    a[[layer]]$start[1] <- 0.5
    expect_error(validate_ichor_sample(a), class = "ichorviz_validation_error")
    a <- example_sample("a")
    a[[layer]]$end[1] <- 3e8
    expect_error(validate_ichor_sample(a), class = "ichorviz_validation_error")
    a <- example_sample("a")
    a[[layer]]$chr[1] <- NA_character_
    expect_error(validate_ichor_sample(a), class = "ichorviz_validation_error")
    a <- example_sample("a")
    a[[layer]]$start[2] <- a[[layer]]$end[1]
    expect_error(validate_ichor_sample(a), class = "ichorviz_validation_error")
  }
  a <- example_sample("a")
  a$bins$logR[] <- NA_real_
  expect_invisible(validate_ichor_sample(a))
  expect_error(plot_ichor_profile(a), "No finite")
  a$bins$copy_number[1] <- -1
  expect_error(validate_ichor_sample(a), class = "ichorviz_validation_error")
})

test_that("component schema ambiguity is an error, not a heuristic", {
  f <- tempfile(fileext = ".cna.seg")
  d <- data.table::fread(fixture_path("example-a.cna.seg"), data.table = FALSE)
  d[["other.logR"]] <- d[["example-a.logR"]]
  data.table::fwrite(d, f, sep = "\t")
  expect_error(read_ichor_cna(f), "Multiple sample prefixes")
  d[["other.logR"]] <- NULL
  d[["logR"]] <- d[["example-a.logR"]]
  data.table::fwrite(d, f, sep = "\t")
  expect_error(read_ichor_cna(f), "Ambiguous")
  d[["logR"]] <- NULL
  names(d) <- sub("example-a.", "", names(d), fixed = TRUE)
  data.table::fwrite(d, f, sep = "\t")
  expect_equal(read_ichor_cna(f)$copy_number, d$copy.number)
})

test_that("seg.txt headers are distinct from canonical seg and normalize", {
  f <- tempfile(fileext = ".seg.txt")
  d <- data.table::fread(fixture_path("example-a.seg"), data.table = FALSE)
  names(d) <- c("ID", "chrom", "start", "end", "call", "copy.number", "num.mark", "seg.median.logR")
  data.table::fwrite(d, f, sep = "\t")
  actual <- read_ichor_segments(f)
  expect_equal(actual$median, read_ichor_segments(fixture_path("example-a.seg"))$median)
  expect_equal(actual$bins, d$num.mark)
})

test_that("parameter parsing rejects conflicts and recognizes candidate diagnostics", {
  f <- tempfile(fileext = ".params.txt")
  lines <- readLines(fixture_path("example-a.params.txt"))
  writeLines(c(lines, "init\tn_est\tphi_est\tBIC", "n0.5-p2\t0.5\t2\tNA"), f)
  expect_equal(read_ichor_params(f)$tumor_fraction, 0.125)
  writeLines(c(lines, "Tumor Fraction:\t0.7"), f)
  expect_error(read_ichor_params(f), "Conflicting")
  writeLines(c(lines, "Tumor Fraction:\tbad"), f)
  expect_error(read_ichor_params(f), "Malformed")
  writeLines(character(), f)
  expect_error(read_ichor_params(f), "No recognized")
  writeLines(c("Sample\tTumor Fraction\tPloidy", "a\t0.1\t2", "b\t0.2\t3"), f)
  expect_error(read_ichor_params(f), "Multiple sample")
  writeLines(c("example-a", "Tumor Fraction:\t0.125", "Ploidy:\t2.1"), f)
  expect_equal(read_ichor_params(f)$sample_id, "example-a")
})

test_that("diagnostic inverse-copy-number infinities are distinct from logR", {
  f <- tempfile(fileext = ".cna.seg")
  d <- data.table::fread(fixture_path("example-a.cna.seg"), data.table = FALSE)
  d[["example-a.logR_Copy_Number"]][1] <- Inf
  data.table::fwrite(d, f, sep = "\t")
  expect_identical(read_ichor_cna(f)$logR_copy_number[1], Inf)
  d[["example-a.logR"]][1] <- Inf
  data.table::fwrite(d, f, sep = "\t")
  expect_error(read_ichor_cna(f), class = "ichorviz_schema_error")
})

test_that("terminal clipping is explicit, recorded, and does not select a build", {
  f <- tempfile(fileext = ".cna.seg")
  d <- data.table::fread(fixture_path("example-a.cna.seg"), data.table = FALSE)[1, ]
  d$start <- 248000001; d$end <- 249000000
  data.table::fwrite(d, f, sep = "\t")
  expect_error(read_ichor_sample(f, genome_build = "hg38", bounds = "error"), "out-of-build")
  expect_warning(a <- read_ichor_sample(f, genome_build = "hg38", bounds = "trim"), "Trimmed 1")
  expect_equal(a$bins$end, 248956422)
  expect_equal(a$coordinate_changes$original_end, 249000000)
  expect_identical(a$coordinate_policy, "trim")
  d$start <- 249000001; d$end <- 250000000
  data.table::fwrite(d, f, sep = "\t")
  expect_error(read_ichor_sample(f, genome_build = "hg38", bounds = "trim"), "out-of-build")
})

test_that("fingerprints are immutable and path retention is opt-in", {
  f <- tempfile(fileext = ".cna.seg")
  file.copy(fixture_path("example-a.cna.seg"), f)
  a <- read_ichor_sample(f, genome_build = "hg38")
  before <- ichor_provenance(a)
  write("\n", file = f, append = TRUE)
  expect_identical(ichor_provenance(a), before)
  expect_false(identical(before$md5, unname(tools::md5sum(f))))
  b <- read_ichor_sample(f, genome_build = "hg38", retain_paths = TRUE)
  expect_equal(ichor_provenance(b)$path, normalizePath(f))
  f2 <- tempfile(); saveRDS(a, f2)
  expect_identical(readRDS(f2), a)
})

test_that("overlap aggregation preserves category semantics and coverage", {
  c <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  a <- c$samples[[1]]
  a$bins <- a$bins[1:2, ]
  a$bins$start <- c(1, 500001); a$bins$end <- c(500000, 1000000)
  a$bins$corrected_call <- c("HETD", "GAIN")
  a$bins$logR <- c(-1, 1)
  c$samples[[1]] <- a
  calls <- ichor_matrix(c, value = "call", chromosomes = "1")
  expect_true(is.na(calls$values[1, 1]))
  expect_true(calls$mixed[1, 1])
  expect_equal(calls$coverage[1, 1], 1)
  continuous <- ichor_matrix(c, value = "logR", chromosomes = "1")
  expect_equal(continuous$values[1, 1], 0)
  c$samples[[1]]$bins$end[1] <- 700000
  c$samples[[1]]$bins$start[2] <- 700001
  calls <- ichor_matrix(c, value = "call", chromosomes = "1")
  expect_equal(calls$values[1, 1], -1)
  expect_true(calls$mixed[1, 1])
  c$samples[[1]]$bins <- a$bins[1, ]
  partial <- ichor_matrix(c, value = "call", chromosomes = "1")
  expect_true(is.na(partial$values[1, 1]))
  expect_equal(partial$coverage[1, 1], 0.5)
  allowed <- ichor_matrix(c, value = "call", chromosomes = "1", min_coverage = 0.5)
  expect_equal(allowed$values[1, 1], -1)
  c$samples[[1]]$bins$corrected_call <- NA_character_
  missing <- ichor_matrix(c, value = "call", chromosomes = "1", min_coverage = 0)
  expect_true(is.na(missing$values[1, 1]))
  expect_equal(missing$coverage[1, 1], 0)
})

test_that("splitting source bins and chromosome-end widths are correct", {
  c <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  m <- ichor_matrix(c, value = "corrected_copy_number", bin_size = 500000, chromosomes = "1")
  expect_equal(unname(m$values[1, 1:8]), rep(c(2, 2, 3, 3), each = 2))
  expect_equal(unname(m$coverage[1, 1:8]), rep(1, 8))
  s <- c$samples[[1]]
  s$bins <- s$bins[1, ]; s$bins$start <- 248000001; s$bins$end <- 248956422
  c$samples[[1]] <- s
  m <- ichor_matrix(c, value = "logR", chromosomes = "1")
  expect_equal(unname(m$coverage[1, ncol(m$values)]), 1)
  expect_error(ichor_matrix(c, value = "call", bin_size = 1), "max_cells")
  expect_error(ichor_matrix(c, value = "call", bin_size = 1.5), "integer")
})

test_that("cohorts fail atomically and enforce metadata order", {
  tab <- data.frame(sample_id = c("001", "002"),
    cna_seg = c(fixture_path("example-a.cna.seg"), fixture_path("example-b.cna.seg")), annotation = c("A", "B"))
  a <- read_ichor_cohort(tab, "hg38", metadata_columns = character())
  expect_identical(names(a$metadata), "sample_id")
  a$metadata <- a$metadata[2:1, , drop = FALSE]
  expect_error(validate_ichor_cohort(a), "order disagree")
  expect_error(read_ichor_cohort(tab, "hg38", workers = 1.5), "integer")
  tab$cna_seg[2] <- "nonexistent.cna.seg"
  expect_error(read_ichor_cohort(tab, "hg38"), "require root")
  expect_error(read_ichor_cohort(tab, "hg38", root = tempdir()), "Cohort import aborted")
  if (.Platform$OS.type != "windows") {
    expect_error(read_ichor_cohort(tab, "hg38", root = tempdir(), workers = 2), "Cohort import aborted")
    c1 <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
    c2 <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38", workers = 2)
    expect_identical(c1, c2)
  }
})

test_that("plotting never changes bin midpoints or mixes optional schemas", {
  a <- example_sample("a"); b <- example_sample("b")
  b$bins$logR_copy_number <- NULL
  b$segments$event <- NULL
  expect_s3_class(ggplot2::ggplot_build(plot_ichor_compare(list(a, b))), "ggplot_built")
  p <- plot_ichor_region(a, "chr1:750000-1500000")
  expect_equal(p$data$x[1], 0.5000005)
  segments <- p$layers[[3]]$data
  expect_gte(min(segments$x), 0.75)
  expect_lte(max(segments$xend), 1.5)
  a$bins$corrected_call <- NULL
  expect_error(plot_ichor_profile(a), "Requested call_column")
  expect_s3_class(plot_ichor_profile(a, call_column = "event"), "ggplot")
  for (bad in c("chr1:", "chr1:1-", "chr1:1.5-2", "chr1:1-2-")) {
    expect_error(parse_ichor_region(bad, "hg38"), class = "ichorviz_region_error")
  }
})

test_that("heatmap drawing succeeds and ordering and missingness are deliberate", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  c <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  m <- ichor_matrix(c, value = "call", chromosomes = "1")
  p <- tempfile(fileext = ".pdf")
  grDevices::pdf(p)
  tryCatch({
    h <- plot_ichor_heatmap(m, annotation_columns = "condition", row_order = rev(m$samples$sample_id),
                            show_row_names = FALSE)
    h <- ComplexHeatmap::draw(h)
    expect_equal(ComplexHeatmap::row_order(h), 2:1)
    m$values[2, ] <- NA_real_
    expect_error(plot_ichor_heatmap(m, cluster_rows = TRUE), "at least two bins")
    expect_s4_class(ComplexHeatmap::draw(plot_ichor_heatmap(m)), "HeatmapList")
    m$samples <- m$samples[2:1, ]
    expect_error(plot_ichor_heatmap(m), "metadata or bin order")
  }, finally = grDevices::dev.off())
  expect_gt(file.info(p)$size, 0)
})
