#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  cat("Usage:\n")
  cat("  Rscript plot_fst_with_outliers.R <windows.tsv> <outliers.tsv> <chromosome> <output.png>\n")
  quit(status = 1)
}

windows_file <- args[1]
outliers_file <- args[2]
chrom_target <- args[3]
output_png <- args[4]

df <- read_tsv(windows_file, col_types = cols()) %>%
  filter(CHROM == chrom_target)

outs <- read_tsv(outliers_file, col_types = cols())

if (nrow(df) == 0) {
  stop("No windows found for chromosome: ", chrom_target)
}

p <- ggplot(df, aes(x = window_start / 1e6, y = fst_mean)) +
  geom_line(color = "grey60") +
  geom_point(color = "grey60", size = 1) +
  geom_point(
    data = outs,
    aes(x = window_start / 1e6, y = fst_mean),
    color = "red",
    size = 2
  ) +
  labs(
    title = paste("FST landrace vs modern -", chrom_target),
    x = "Position (Mb)",
    y = "Mean FST per window"
  ) +
  theme_bw(base_size = 14)

ggsave(output_png, plot = p, width = 10, height = 4, dpi = 300)

cat("Outlier plot written to: ", output_png, "\n", sep = "")
