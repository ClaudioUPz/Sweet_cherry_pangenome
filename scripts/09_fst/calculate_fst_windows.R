#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage:\n")
  cat("  Rscript calculate_fst_windows.R \\\n")
  cat("    <fst_table.tsv> \\\n")
  cat("    <window_size_bp> \\\n")
  cat("    <output_windows.tsv>\n\n")
  quit(status = 1)
}

in_fst <- args[1]
window_size <- as.numeric(args[2])
out_file <- args[3]

cat("Input FST file: ", in_fst, "\n", sep = "")
cat("Window size (bp): ", window_size, "\n", sep = "")
cat("Output file: ", out_file, "\n\n", sep = "")

fst <- read_tsv(in_fst, col_types = cols())

cat("Columns detected in FST file:\n")
print(names(fst))

fst_col_candidates <- c("Fst_Hudson", "Fst", "fst", "FST")
fst_col <- fst_col_candidates[fst_col_candidates %in% names(fst)][1]

if (is.na(fst_col)) {
  stop("No FST column found. Expected one of: Fst_Hudson, Fst, fst, FST.")
}

if (!("CHROM" %in% names(fst)) || !("POS" %in% names(fst))) {
  stop("Required columns CHROM and POS were not found in the FST table.")
}

cat("\nUsing FST column: ", fst_col, "\n\n", sep = "")

fst_clean <- fst %>%
  filter(!is.na(.data[[fst_col]]))

cat("Total rows in FST file: ", nrow(fst), "\n", sep = "")
cat("Rows with non-missing FST values: ", nrow(fst_clean), "\n\n", sep = "")

fst_windows <- fst_clean %>%
  mutate(
    window_start = floor((POS - 1) / window_size) * window_size + 1,
    window_end   = window_start + window_size - 1
  ) %>%
  group_by(CHROM, window_start, window_end) %>%
  summarise(
    n_snps     = n(),
    fst_mean   = mean(.data[[fst_col]], na.rm = TRUE),
    fst_median = median(.data[[fst_col]], na.rm = TRUE),
    fst_max    = max(.data[[fst_col]], na.rm = TRUE),
    fst_min    = min(.data[[fst_col]], na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(CHROM, window_start)

cat("Number of windows generated: ", nrow(fst_windows), "\n", sep = "")

write_tsv(fst_windows, out_file)

cat("\nFST windows written to: ", out_file, "\n", sep = "")
