#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage:\n")
  cat("  Rscript get_fst_outliers_all_chr.R <windows.tsv> <output_prefix>\n")
  quit(status = 1)
}

windows_file <- args[1]
out_prefix <- args[2]

df_all <- read_tsv(windows_file, col_types = cols())

chroms <- sort(unique(df_all$CHROM))

cat("Detected chromosomes: ", paste(chroms, collapse = ", "), "\n", sep = "")

for (chrom_target in chroms) {
  df <- df_all %>%
    filter(CHROM == chrom_target)

  if (nrow(df) == 0) next

  threshold <- quantile(df$fst_mean, 0.95, na.rm = TRUE)

  cat("\n", chrom_target, " - FST threshold (95th percentile): ", threshold, "\n", sep = "")

  outliers <- df %>%
    filter(fst_mean >= threshold)

  output_file <- paste0(out_prefix, "_", chrom_target, ".tsv")

  write_tsv(outliers, output_file)

  cat("Outliers written to: ", output_file,
      " | n_outliers = ", nrow(outliers), "\n", sep = "")
}
