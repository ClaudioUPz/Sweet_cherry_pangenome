#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage:\n",
    "Rscript 06_plot_go_chr_dominance.R ",
    "<go_by_chr.tsv> <go_enrichment_named.tsv> <output_prefix>\n"
  )
}

bychr_file <- args[1]
go_file <- args[2]
output_prefix <- args[3]

cat("Reading input files...\n")

bychr <- read_tsv(bychr_file, show_col_types = FALSE)
go <- read_tsv(go_file, show_col_types = FALSE)

required_bychr <- c("GO", "CHR", "N_genes")
missing_bychr <- setdiff(required_bychr, colnames(bychr))
if (length(missing_bychr) > 0) {
  stop("Missing required columns in chromosome summary file: ", paste(missing_bychr, collapse = ", "))
}

required_go <- c("ID", "Description", "p.adjust", "Count")
missing_go <- setdiff(required_go, colnames(go))
if (length(missing_go) > 0) {
  stop("Missing required columns in GO enrichment file: ", paste(missing_go, collapse = ", "))
}

df <- bychr %>%
  left_join(
    go %>% select(ID, Description, p.adjust, Count),
    by = c("GO" = "ID")
  ) %>%
  mutate(
    Description = ifelse(is.na(Description) | Description == "", GO, Description)
  ) %>%
  group_by(GO) %>%
  mutate(
    total_genes = sum(N_genes),
    prop = N_genes / total_genes
  ) %>%
  ungroup() %>%
  mutate(
    GO_label = paste0(Description, " (n=", total_genes, ")")
  )

go_order <- df %>%
  distinct(GO, GO_label, p.adjust) %>%
  arrange(p.adjust) %>%
  pull(GO_label)

df$GO_label <- factor(df$GO_label, levels = rev(unique(go_order)))

chr_info <- df %>%
  distinct(CHR) %>%
  mutate(
    chr_num = suppressWarnings(as.numeric(str_extract(CHR, "\\d+")))
  ) %>%
  arrange(chr_num, CHR)

df$CHR <- factor(df$CHR, levels = chr_info$CHR)

cat("Generating stacked proportion plot...\n")

p <- ggplot(df, aes(x = prop, y = GO_label, fill = CHR)) +
  geom_col(width = 0.75) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Proportion of genes",
    y = NULL,
    fill = "Chromosome"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

ggsave(
  paste0(output_prefix, ".stacked_prop_named.png"),
  p,
  width = 9,
  height = 4.8,
  dpi = 300
)

ggsave(
  paste0(output_prefix, ".stacked_prop_named.pdf"),
  p,
  width = 9,
  height = 4.8
)

cat("Done\n")
cat("Saved:\n")
cat(" - ", paste0(output_prefix, ".stacked_prop_named.png"), "\n", sep = "")
cat(" - ", paste0(output_prefix, ".stacked_prop_named.pdf"), "\n", sep = "")
