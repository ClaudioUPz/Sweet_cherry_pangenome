#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage:\n",
    "Rscript 07_go_windows_positions.R ",
    "<go_genes_in_outlier_windows.tsv> <output_prefix>\n"
  )
}

input_file <- args[1]
output_prefix <- args[2]

df <- read_tsv(input_file, col_names = FALSE, show_col_types = FALSE)

if (ncol(df) != 8) {
  stop("Input file must contain exactly 8 columns.")
}

colnames(df) <- c(
  "gene_chr",
  "gene_start",
  "gene_end",
  "gene",
  "GO",
  "window_chr",
  "window_start",
  "window_end"
)

chr_mismatch <- df %>%
  filter(gene_chr != window_chr)

if (nrow(chr_mismatch) > 0) {
  warning("Some rows contain different gene and window chromosome values. Using window chromosome for summarization.")
}

windows <- df %>%
  mutate(chr = window_chr) %>%
  group_by(GO, chr, window_start, window_end) %>%
  summarise(
    n_genes = n(),
    genes = paste(unique(gene), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(GO, chr, window_start)

write_tsv(windows, paste0(output_prefix, ".go_windows.tsv"))

p <- ggplot(
  windows,
  aes(
    x = window_start / 1e6,
    xend = window_end / 1e6,
    y = chr,
    yend = chr
  )
) +
  geom_segment(linewidth = 3, color = "red") +
  facet_wrap(~GO, scales = "free_x") +
  labs(
    x = "Position (Mb)",
    y = "Chromosome",
    title = "FST outlier windows contributing genes to enriched GO terms"
  ) +
  theme_bw()

ggsave(
  paste0(output_prefix, ".go_windows_plot.png"),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  paste0(output_prefix, ".go_windows_plot.pdf"),
  p,
  width = 10,
  height = 6
)

cat("Saved:\n")
cat(" - ", paste0(output_prefix, ".go_windows.tsv"), "\n", sep = "")
cat(" - ", paste0(output_prefix, ".go_windows_plot.png"), "\n", sep = "")
cat(" - ", paste0(output_prefix, ".go_windows_plot.pdf"), "\n", sep = "")
