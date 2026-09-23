# Synthetic documentation data, not patient data or fitted ichorCNA output.
# Source this file, then call write_demo_cohort() to create temporary inputs.
write_demo_cohort <- function(directory = tempfile("seena-demo-")) {
  if (file.exists(directory) && (!dir.exists(directory) ||
      length(list.files(directory, all.files = TRUE, no.. = TRUE)))) {
    stop("Choose a new or empty directory for the synthetic demo.")
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(directory)) stop("Cannot create the demo directory.")

  # Generating documentation must not change the caller's random-number stream.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(41027, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  layout <- seeNA::ichor_genome_layout("hg38", c(as.character(1:22), "X"))
  grid <- do.call(rbind, lapply(seq_len(nrow(layout)), function(k) {
    start <- seq(1, layout$length[k], by = 5e5)
    data.frame(chr = layout$chr[k], start = start,
               end = pmin(start + 5e5 - 1, layout$length[k]))
  }))
  chr_index <- match(grid$chr, layout$chr)
  mid <- (grid$start + grid$end) / 2
  fraction <- mid / layout$length[chr_index]
  templates <- rbind(
    A = c(3, 2, 1, 2, 3, 2, 3, 4, 2, 1, 2, 3, 1, 2, 2, 1, 2, 1, 2, 3, 2, 2, 3),
    B = c(2, 3, 2, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2, 3, 1, 2, 3, 2, 1, 2, 3, 2, 2),
    C = c(1, 2, 3, 3, 2, 1, 4, 2, 3, 2, 1, 2, 3, 2, 3, 2, 2, 3, 2, 1, 2, 3, 1)
  )
  pair_tf <- c(0.46, 0.32, 0.58, 0.20, 0.50, 0.38, 0.64, 0.26, 0.55, 0.42, 0.30, 0.60)
  tf <- round(as.vector(rbind(pair_tf, pair_tf * 0.6 + 0.02)), 3)
  manifest <- data.frame(sample_id = sprintf("demo-%02d", 1:24),
    group = rep(c("A", "B", "C"), each = 8),
    pair = rep(sprintf("pair-%02d", 1:12), each = 2),
    timepoint = rep(c("T1", "T2"), 12), stringsAsFactors = FALSE)
  manifest[["TF (%)"]] <- tf * 100
  manifest$cna_seg <- paste0(manifest$sample_id, ".cna.seg")
  manifest$seg <- paste0(manifest$sample_id, ".seg")
  manifest$params <- paste0(manifest$sample_id, ".params.txt")
  write_tsv <- function(x, name) utils::write.table(x, file.path(directory, name),
    sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

  for (i in seq_len(nrow(manifest))) {
    pair <- ceiling(i / 2)
    copy <- rep(2, nrow(grid))
    lo <- 0.08 + ((chr_index * 3 + pair) %% 5) * 0.045
    hi <- 0.58 + ((chr_index + pair) %% 4) * 0.08
    altered <- fraction >= lo & fraction < hi
    copy[altered] <- templates[manifest$group[i], chr_index[altered]]
    focal <- (chr_index + pair) %% 6 == 0 & fraction >= 0.82 & fraction < 0.9
    copy[focal] <- c(1, 3, 4)[(chr_index[focal] + pair) %% 3 + 1]
    if (pair %% 3 == 1) copy[grid$chr == "17" & mid >= 20e6 & mid < 28e6] <- 0
    if (i %% 2 == 0) copy[chr_index == pair + 2 & fraction > 0.4 & fraction < 0.65] <- 2

    # The first pair deliberately exercises each directional agreement category.
    if (pair == 1) {
      copy[grid$chr == "1"] <- 2
      copy[grid$chr == "1" & mid >= 25e6 & mid < 55e6] <- 3
      copy[grid$chr == "1" & mid >= 75e6 & mid < 95e6] <- if (i == 1) 1 else 2
      copy[grid$chr == "1" & mid >= 115e6 & mid < 140e6] <- if (i == 1) 2 else 4
      copy[grid$chr == "1" & mid >= 165e6 & mid < 190e6] <- if (i == 1) 3 else 1
    }
    keep <- !(grid$chr == "9" & fraction > 0.44 & fraction < 0.49)
    keep <- keep & !(chr_index == (i * 3) %% 22 + 1 & fraction > 0.7 & fraction < 0.78)
    if (pair == 1) keep <- keep & !(grid$chr == "1" & mid >= 202e6 & mid < 210e6)
    d <- grid[keep, ]
    cn <- copy[keep]
    call <- c("HOMD", "HETD", "NEUT", "GAIN", "AMP")[cn + 1]
    logr <- round(log2((2 * (1 - tf[i]) + cn * tf[i]) / 2) +
                    stats::rnorm(nrow(d), sd = 0.055), 6)
    id <- manifest$sample_id[i]
    bins <- data.frame(d, copy.number = cn, event = call, logR = logr,
                       Corrected_Copy_Number = cn, Corrected_Call = call)
    names(bins)[4:8] <- paste(id, names(bins)[4:8], sep = ".")
    write_tsv(bins, manifest$cna_seg[i])

    # Split at gaps as well as state boundaries; missing intervals stay absent.
    key <- paste(d$chr, cn)
    new_segment <- c(TRUE, key[-1] != key[-length(key)] | d$start[-1] != utils::head(d$end, -1) + 1)
    rows <- split(seq_len(nrow(d)), cumsum(new_segment))
    first <- vapply(rows, function(j) j[1], integer(1))
    last <- vapply(rows, function(j) j[length(j)], integer(1))
    segments <- data.frame(sample = id, chr = d$chr[first], start = d$start[first],
      end = d$end[last], event = call[first], copy.number = cn[first], bins = lengths(rows),
      median = vapply(rows, function(j) round(stats::median(logr[j]), 6), numeric(1)))
    write_tsv(segments, manifest$seg[i])
    params <- data.frame(Sample = id, "Tumor Fraction" = tf[i], Ploidy = 2, check.names = FALSE)
    write_tsv(params, manifest$params[i])
  }
  path <- file.path(directory, "manifest.csv")
  utils::write.csv(manifest, path, row.names = FALSE)
  path
}
