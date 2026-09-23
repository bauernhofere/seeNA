test_that("the package namespace, fixtures and conditions use the new name", {
  expect_identical(unname(getNamespaceName(environment(ichor_matrix))), "seeNA")
  expect_true(file.exists(system.file("extdata", "example-a.seg", package = "seeNA")))
  err <- tryCatch(.ichor_abort("test-only error"), error = identity)
  expect_s3_class(err, "seena_error")
})
