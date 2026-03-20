#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage:\n")
  cat("  Rscript plot_fst_windows.R \\\n")
  cat("    <windows.tsv> \\\n")
  cat("    <chromosome> \\\n")
  cat("    <output.png>\n\n")
  quit(status = 1)
}

windows_file <- args[1]
chrom_target <- args[2]
output_png <- args[3]

cat("Input windows file: ", windows_file, "\n", sep = "")
cat("Target chromosome: ", chrom_target, "\n", sep = "")
cat("Output plot: ", output_png, "\n\n", sep = "")

df <- read_tsv(windows_file, col_types = cols())

df_chr <- df %>%
  filter(CHROM == chrom_target)

if (nrow(df_chr) == 0) {
  stop("No windows found for chromosome: ", chrom_target)
}

p <- ggplot(df_chr, aes(x = window_start / 1e6, y = fst_mean)) +
  geom_line(color = "steelblue", linewidth = 0.7) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    title = paste("FST landrace vs modern -", chrom_target),
    x = "Position (Mb)",
    y = "Mean FST per window"
  ) +
  theme_bw(base_size = 14)

ggsave(output_png, plot = p, width = 10, height = 4, dpi = 300)

cat("Plot written to: ", output_png, "\n", sep = "")
