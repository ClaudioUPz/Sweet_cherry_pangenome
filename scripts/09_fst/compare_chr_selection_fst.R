#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage:\n")
  cat("Rscript compare_chr_selection_fst.R <windows.tsv> <output_prefix>\n")
  quit(status = 1)
}

windows_file <- args[1]
output_prefix <- args[2]

df <- read_tsv(windows_file, col_types = cols())

required_cols <- c("CHROM", "fst_mean")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
}

threshold <- quantile(df$fst_mean, 0.99, na.rm = TRUE)
cat("Genome-wide FST threshold (99th percentile):", threshold, "\n")

summary_chr <- df %>%
  group_by(CHROM) %>%
  summarise(
    total_windows = n(),
    n_top = sum(fst_mean >= threshold, na.rm = TRUE),
    prop_top = n_top / total_windows,
    mean_fst = mean(fst_mean, na.rm = TRUE),
    max_fst = max(fst_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    CHROM_clean = sub("__.*", "", CHROM),
    chrom_num = suppressWarnings(as.numeric(sub("^chr", "", CHROM_clean)))
  ) %>%
  arrange(chrom_num, CHROM_clean)

summary_chr$CHROM_clean <- factor(summary_chr$CHROM_clean, levels = unique(summary_chr$CHROM_clean))

write_tsv(summary_chr, paste0(output_prefix, "_chr_summary_p99.tsv"))
cat("Saved summary table\n")
print(summary_chr %>% arrange(desc(n_top)))

p1 <- ggplot(summary_chr, aes(x = CHROM_clean, y = n_top)) +
  geom_col() +
  labs(
    title = "Number of high-FST windows per chromosome (p99)",
    subtitle = paste0("Genome-wide threshold = ", round(threshold, 4)),
    x = "Chromosome",
    y = "Number of windows (FST >= p99)"
  ) +
  theme_bw(base_size = 14)

ggsave(
  paste0(output_prefix, "_nTop_barplot.png"),
  plot = p1,
  width = 8,
  height = 4,
  dpi = 300
)

p2 <- ggplot(summary_chr, aes(x = CHROM_clean, y = prop_top)) +
  geom_col() +
  labs(
    title = "Proportion of high-FST windows per chromosome (p99)",
    x = "Chromosome",
    y = "Proportion of windows (FST >= p99)"
  ) +
  theme_bw(base_size = 14)

ggsave(
  paste0(output_prefix, "_propTop_barplot.png"),
  plot = p2,
  width = 8,
  height = 4,
  dpi = 300
)

p3 <- ggplot(summary_chr, aes(x = CHROM_clean, y = max_fst)) +
  geom_col() +
  labs(
    title = "Maximum FST per chromosome",
    x = "Chromosome",
    y = "Maximum FST"
  ) +
  theme_bw(base_size = 14)

ggsave(
  paste0(output_prefix, "_maxFst_barplot.png"),
  plot = p3,
  width = 8,
  height = 4,
  dpi = 300
)

cat("Saved plots\n")
