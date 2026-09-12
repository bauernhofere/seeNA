# Genome coordinates and regions ------------------------------------------

.chromosome_sizes <- function(build) {
  build <- match.arg(build, c("hg19", "hg38"))
  sizes <- if (build == "hg19") {
    c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
      159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
      115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
      59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  } else {
    c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
      159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
      114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
      58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  }
  data.frame(chr = .chr_levels, length = sizes, stringsAsFactors = FALSE)
}

#' Reference chromosome lengths and cumulative offsets
#' @param genome_build Explicit hg19 or hg38.
#' @param chromosomes Chromosomes to include; default autosomes and X.
#'   Selected chromosomes are always ordered genomically. Include Y explicitly.
#' @return Data frame with chr, length, offset, mid and boundary in base pairs.
#'   Offsets depend on the chromosome set; reuse the same layout across samples.
#' @export
ichor_genome_layout <- function(genome_build, chromosomes = c(as.character(1:22), "X")) {
  d <- .chromosome_sizes(genome_build)
  chromosomes <- unique(.normalize_chr(chromosomes))
  if (!length(chromosomes) || anyNA(chromosomes) || !all(chromosomes %in% d$chr)) .ichor_abort("Unsupported chromosomes.")
  d <- d[d$chr %in% chromosomes, , drop = FALSE]
  d$offset <- c(0, utils::head(cumsum(d$length), -1))
  d$mid <- d$offset + d$length / 2
  d$boundary <- d$offset + d$length
  d
}

.add_genome_coordinates <- function(d, layout, segments = FALSE) {
  i <- match(d$chr, layout$chr)
  if (anyNA(i)) .ichor_abort("Data contain chromosomes absent from the genome layout.")
  if (segments) {
    d$x <- layout$offset[i] + d$start
    d$xend <- layout$offset[i] + d$end
  } else {
    d$x <- layout$offset[i] + (d$start + d$end) / 2
  }
  d
}

#' Parse a genomic region
#'
#' Accepts a whole chromosome such as `"chr8"`, or an interval such as
#' `"chr8:117000000-138000000"`. Commas in coordinates are allowed.
#'
#' @param region Region string.
#' @param genome_build Genome build used for bounds checking.
#' @return A one-row data frame with `chr`, `start`, and `end`.
#' @export
parse_ichor_region <- function(region, genome_build) {
  genome_build <- match.arg(genome_build, c("hg19", "hg38"))
  if (length(region) != 1L || is.na(region) || !nzchar(region)) {
    .ichor_abort("region must be a non-empty scalar.", "ichorviz_region_error")
  }
  clean <- gsub("[,]", "", trimws(region))
  if (!grepl("^(chr)?([0-9]+|X|Y)(:[0-9]+-[0-9]+)?$", clean, ignore.case = TRUE)) {
    .ichor_abort("Use chromosome or chr:start-end with integer coordinates.", "ichorviz_region_error")
  }
  pieces <- strsplit(clean, ":", fixed = TRUE)[[1]]
  chr <- .normalize_chr(pieces[1])
  sizes <- .chromosome_sizes(genome_build)
  chr_length <- sizes$length[match(chr, sizes$chr)]
  if (is.na(chr_length)) .ichor_abort(sprintf("Unsupported chromosome: %s", chr), "ichorviz_region_error")

  if (length(pieces) == 1L) return(data.frame(chr = chr, start = 1, end = chr_length))
  if (length(pieces) != 2L) .ichor_abort("Use region syntax chr:start-end.", "ichorviz_region_error")
  bounds <- strsplit(pieces[2], "-", fixed = TRUE)[[1]]
  if (length(bounds) != 2L) .ichor_abort("Use region syntax chr:start-end.", "ichorviz_region_error")
  start <- .as_number(bounds[1]); end <- .as_number(bounds[2])
  if (!is.finite(start) || !is.finite(end) || start < 1 || end < start || end > chr_length) {
    .ichor_abort(sprintf("Region lies outside %s (length %s).", chr,
                         format(chr_length, big.mark = ",", scientific = FALSE)),
                 "ichorviz_region_error")
  }
  data.frame(chr = chr, start = start, end = end)
}

.subset_interval <- function(d, region) {
  keep <- d$chr == region$chr & d$end >= region$start & d$start <= region$end
  out <- d[keep, , drop = FALSE]
  if (nrow(out)) {
    out$start <- pmax(out$start, region$start)
    out$end <- pmin(out$end, region$end)
  }
  out
}
