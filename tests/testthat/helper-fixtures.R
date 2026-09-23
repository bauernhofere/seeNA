fixture_path <- function(name) {
  system.file("extdata", name, package = "seeNA", mustWork = TRUE)
}

example_sample <- function(id = c("a", "b")) {
  id <- match.arg(id)
  read_ichor_sample(
    fixture_path(sprintf("example-%s.cna.seg", id)),
    fixture_path(sprintf("example-%s.seg", id)),
    fixture_path(sprintf("example-%s.params.txt", id)),
    genome_build = "hg38"
  )
}
