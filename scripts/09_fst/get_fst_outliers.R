#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage:\n")
  cat("  Rscript get_fst_outliers.R <windows.tsv> <chromosome> <output.tsv>\n")
  quit(status = 1)
}

windows_file <- args[1]
chrom_target <- args[2]
output_file <- args[3]

df <- read_tsv(windows_file, col_types = cols()) %>%
  filter(CHROM == chrom_target)

if (nrow(df) == 0) {
  stop("No windows found for chromosome: ", chrom_target)
}

threshold <- quantile(df$fst_mean, 0.95, na.rm = TRUE)

cat("FST threshold (95th percentile): ", threshold, "\n", sep = "")

outliers <- df %>%
  filter(fst_mean >= threshold)

write_tsv(outliers, output_file)

cat("Outliers written to: ", output_file, "\n", sep = "")
