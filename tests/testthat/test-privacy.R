test_that("public fixtures contain no project patient identifiers", {
  files <- list.files(system.file("extdata", package = "seeNA"), full.names = TRUE)
  text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")

  expect_false(grepl("P[0-9]{3}", text))
  expect_false(grepl("/Users/", text, fixed = TRUE))
})
