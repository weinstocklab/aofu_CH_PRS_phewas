#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: build_phecode_resources.R <expanded_map.tsv> <phecode_info.tsv>")
}

expanded_out <- args[[1]]
info_out <- args[[2]]

download_rda <- function(name) {
  url <- sprintf("https://raw.githubusercontent.com/PheWAS/PheWAS/master/data/%s.rda", name)
  dest <- tempfile(fileext = ".rda")
  download.file(url, destfile = dest, quiet = TRUE)
  loaded <- load(dest)
  get(loaded[[1]])
}

phecode_map <- download_rda("phecode_map")
phecode_rollup_map <- download_rda("phecode_rollup_map")
pheinfo <- download_rda("pheinfo")

expanded <- merge(
  phecode_map,
  phecode_rollup_map,
  by.x = "phecode",
  by.y = "code",
  all = FALSE
)
expanded <- unique(expanded[, c("vocabulary_id", "code", "phecode_unrolled")])
names(expanded) <- c("vocabulary_id", "concept_code", "phecode")

dir.create(dirname(expanded_out), showWarnings = FALSE, recursive = TRUE)
write.table(expanded, expanded_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pheinfo, info_out, sep = "\t", quote = FALSE, row.names = FALSE)
