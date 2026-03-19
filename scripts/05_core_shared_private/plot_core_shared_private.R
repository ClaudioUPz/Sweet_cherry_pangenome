#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

tsv_file <- "per_sample_core_shared_private_v3.tsv"

if (!file.exists(tsv_file)) {
  stop(paste("ERROR: File not found:", tsv_file))
}

df <- read_tsv(tsv_file, col_types = cols())

# --------------------------------
# Rename samples
# --------------------------------
df <- df %>%
  mutate(
    sample = case_when(
      sample == "Regina"    ~ "Regina_hap1",
      sample == "H2Regina"  ~ "Regina_hap2",
      TRUE ~ sample
    ),
    sample = gsub("_h1$", "_hap1", sample),
    sample = gsub("_h2$", "_hap2", sample)
  )

# --------------------------------
# Fractions
# --------------------------------
df2 <- df %>%
  mutate(
    core_frac    = core_bp    / total_bp,
    shared_frac  = shared_bp  / total_bp,
    private_frac = private_bp / total_bp
  )

# --------------------------------
# Long format
# --------------------------------
df_long <- df2 %>%
  select(sample, core_frac, shared_frac, private_frac) %>%
  pivot_longer(
    cols = c(core_frac, shared_frac, private_frac),
    names_to = "category",
    values_to = "fraction"
  ) %>%
  mutate(
    category = recode(category,
                      core_frac    = "Core",
                      shared_frac  = "Shared",
                      private_frac = "Private"),
    category = factor(category, levels = c("Core", "Shared", "Private"))
  )

df_long$sample <- factor(df_long$sample, levels = unique(df_long$sample))

# --------------------------------
# Plot
# --------------------------------
p <- ggplot(df_long, aes(x = sample, y = fraction, fill = category)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Core"="#1b9e77", "Shared"="#7570b3", "Private"="#d95f02")) +
  labs(
    x = "Genome / haplotype",
    y = "Fraction of genomic bases",
    fill = "Category",
    title = "Distribution of core, shared, and private regions per genome"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("Fig1C_core_shared_private.pdf", p, width = 9, height = 7)
ggsave("Fig1C_core_shared_private.png", p, width = 9, height = 7, dpi = 300)

cat("✔ Generated:\n")
cat("   - Fig1C_core_shared_private.pdf\n")
cat("   - Fig1C_core_shared_private.png\n")