#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage:\n")
  cat("Rscript plot_fst_with_outliers_all_chr.R <windows.tsv> <outliers_dir> <output_prefix>\n\n")
  cat("Example:\n")
  cat("Rscript plot_fst_with_outliers_all_chr.R snps.windows.tsv Fst_outliers Fst\n")
  quit(status = 1)
}

windows_file <- args[1]
outliers_dir <- args[2]
output_prefix <- args[3]

outliers_dir <- sub("/+$", "", outliers_dir)

df_all <- read_tsv(windows_file, col_types = cols())

required_cols <- c("CHROM", "window_start", "fst_mean")
missing_cols <- setdiff(required_cols, colnames(df_all))

if (length(missing_cols) > 0) {
  stop("Missing required columns in windows file: ", paste(missing_cols, collapse = ", "))
}

chromosomes <- sort(unique(df_all$CHROM))
cat("Detected chromosomes:", paste(chromosomes, collapse = ", "), "\n")

for (chrom_target in chromosomes) {
  df <- df_all %>% filter(CHROM == chrom_target)

  if (nrow(df) == 0) {
    next
  }

  outlier_file <- file.path(outliers_dir, paste0(output_prefix, "_", chrom_target, ".tsv"))

  if (!file.exists(outlier_file)) {
    cat("Warning: outlier file not found for", chrom_target, "->", outlier_file, "\n")
    next
  }

  outs <- read_tsv(outlier_file, col_types = cols())

  missing_out_cols <- setdiff(required_cols, colnames(outs))
  if (length(missing_out_cols) > 0) {
    stop("Missing required columns in outlier file ", outlier_file, ": ",
         paste(missing_out_cols, collapse = ", "))
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
      title = paste("FST (landrace vs modern breeding) -", chrom_target),
      x = "Position (Mb)",
      y = "Mean FST per 100 kb window"
    ) +
    theme_bw(base_size = 14)

  output_png <- paste0(output_prefix, "_", chrom_target, "_fst_outliers.png")
  ggsave(output_png, plot = p, width = 10, height = 4, dpi = 300)

  cat("Saved:", output_png, "\n")
}
