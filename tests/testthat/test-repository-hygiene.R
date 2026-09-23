test_that("repository ignores private inputs and keeps only explicit public fixtures", {
  root <- testthat::test_path("..", "..")
  skip_if_not(file.exists(file.path(root, ".gitignore")) && file.exists(file.path(root, ".git")),
              "Repository-only check; Git metadata is not shipped in the R package")
  skip_if(!nzchar(Sys.which("git")), "Git is required for the repository-only check")
  ignored <- c(
    ".Renviron", ".env", ".env.local", ".claude/settings.local.json", ".pi/session.json",
    "private/cohort.csv", "data/cohort.csv", "input.xlsx", "input.xls", "manifest.tsv",
    "fitted.seg.txt", "sample.cna.seg", "sample.params.txt", "sample.seg",
    "sample.bam", "sample.cram", "sample.vcf.gz", "cohort.rds", "cohort.RData",
    "output/figure.html", "figure.pdf", "figure.png", "figure.svg", "figure.tiff",
    "seeNA.Rcheck/00check.log", "seeNA_0.0.0.9003.tar.gz", "private-key.pem",
    "inst/extdata/unapproved.cna.seg", "inst/extdata/unapproved.params.txt",
    "inst/extdata/unapproved.seg", "inst/extdata/unapproved.csv"
  )
  public <- c(
    "AGENTS.md", "README.md", "README.Rmd", "DESCRIPTION", "R/seeNA-package.R",
    "inst/extdata/example-a.cna.seg", "inst/extdata/example-a.params.txt", "inst/extdata/example-a.seg",
    "inst/extdata/example-b.cna.seg", "inst/extdata/example-b.params.txt", "inst/extdata/example-b.seg",
    "inst/extdata/example-manifest.csv", "inst/extdata/README.md",
    paste0("docs/figures/README-", c("region", "profile", "heatmap", "concordance"), "-1.png"),
    "tests/testthat/_snaps/snapshots/profile.svg"
  )
  # Disable personal excludes so another contributor gets the same safeguards.
  result <- system2(Sys.which("git"), c("-C", shQuote(normalizePath(root)), "-c",
    shQuote(paste0("core.excludesFile=", tempfile())), "check-ignore", "--no-index", "--stdin"),
    input = c(ignored, public), stdout = TRUE, stderr = TRUE)
  expect_null(attr(result, "status"))
  expect_setequal(result, ignored)
})
