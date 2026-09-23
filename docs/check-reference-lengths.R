# Optional reference audit, not a network-dependent package test.
pkgload::load_all(".", quiet = TRUE)
for (build in c("hg19", "hg38")) {
  url <- sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes", build, build)
  path <- tempfile()
  utils::download.file(url, path, mode = "wb", quiet = TRUE)
  reference <- utils::read.table(path, col.names = c("chr", "length"),
                                 colClasses = c("character", "numeric"))
  layout <- seeNA::ichor_genome_layout(build, c(as.character(1:22), "X", "Y"))
  expected <- reference$length[match(paste0("chr", layout$chr), reference$chr)]
  stopifnot(!anyNA(expected), identical(layout$length, expected))
  cat(sprintf("%s: all 24 primary chromosome lengths match; downloaded file MD5 %s\n",
              build, unname(tools::md5sum(path))))
  unlink(path)
}
