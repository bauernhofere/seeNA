test_that("whole chromosomes and intervals parse", {
  chr <- parse_ichor_region("chr8", "hg38")
  interval <- parse_ichor_region("chr8:1,000,000-2,000,000", "hg38")

  expect_equal(chr$chr, "8")
  expect_equal(chr$end, 145138636)
  expect_equal(interval, data.frame(chr = "8", start = 1000000, end = 2000000))
})

test_that("invalid intervals are rejected", {
  expect_error(parse_ichor_region("chr1:20-10", "hg38"), class = "ichorviz_region_error")
  expect_error(parse_ichor_region("chr27", "hg38"), class = "ichorviz_region_error")
})
