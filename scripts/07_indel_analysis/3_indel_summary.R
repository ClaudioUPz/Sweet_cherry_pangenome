#!/usr/bin/env Rscript

# Read 7.INDEL/chrN/indels_details_chrN.csv and generate:
# - Raw counts of ALT INDELs per sample and chromosome
# - Normalization by FASTA length (Mb) of each sample in that chromosome
# - Per-chromosome and genome-wide plots (robust PNG device)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(viridis)
  library(tidyr)
})
options(dplyr.summarise.inform = FALSE)

BASE <- "/path/to/project"
chr_ids <- paste0("chr", 1:8)   # For testing only chr1: chr_ids <- "chr1"
details_name <- "indels_details_%s.csv"

plots_dir <- file.path(BASE, "7.INDEL", "plots_per_sample")
sums_dir  <- file.path(BASE, "7.INDEL", "summaries_per_sample")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sums_dir, showWarnings = FALSE, recursive = TRUE)

save_png <- function(plot, file, width = 12, height = 8, dpi = 300, bg = "white") {
  png(filename = file, width = width * dpi, height = height * dpi, res = dpi, bg = bg)
  print(plot)
  dev.off()
  message("PNG saved: ", file)
}

# Total FASTA length (bp), ignoring header lines
fasta_bp <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  sum(nchar(lines[!startsWith(lines, ">")]))
}

global_counts <- tibble(SAMPLE = character(), ALT_sites = double())
global_bp     <- tibble(SAMPLE = character(), bp = double())

for (chr in chr_ids) {
  message("\n========== ", chr, " ==========")
  csv <- file.path(BASE, "7.INDEL", chr, sprintf(details_name, chr))
  smp <- file.path(BASE, "7.INDEL", chr, "samples.txt")  # Generated in Step 1

  if (!file.exists(csv)) {
    message("Skipping (file not found): ", csv)
    next
  }

  d <- read_csv(
    csv,
    show_col_types = FALSE,
    col_types = cols(
      CHROM  = col_character(),
      POS    = col_double(),
      REF    = col_character(),
      ALT    = col_character(),
      SAMPLE = col_character(),
      GT     = col_character()
    )
  )

  message("ALT rows in ", chr, ": ", format(nrow(d), big.mark = ","))

  samples_chr <- if (file.exists(smp)) readLines(smp) else sort(unique(d$SAMPLE))

  # Raw count per sample (ALT observations)
  sum_raw <- d %>%
    count(SAMPLE, name = "ALT_sites") %>%
    arrange(desc(ALT_sites))

  # Sequence length per sample in this chromosome
  # FASTA files are expected at: 5.pangenome/chrN/<SAMPLE>_chrN.fa
  fasta_paths <- file.path(BASE, "5.pangenome", chr, paste0(samples_chr, "_", chr, ".fa"))
  bp_vec <- purrr::map_dbl(fasta_paths, fasta_bp)
  len_df <- tibble(SAMPLE = samples_chr, bp = bp_vec)

  # Normalize by Mb
  sum_norm <- sum_raw %>%
    left_join(len_df, by = "SAMPLE") %>%
    mutate(ALT_per_Mb = if_else(!is.na(bp) & bp > 0, ALT_sites / (bp / 1e6), NA_real_))

  # Save summary table for this chromosome
  out_tab <- file.path(sums_dir, paste0("indels_summary_", chr, "_per_sample.csv"))
  write_csv(sum_norm, out_tab)
  message("Table saved: ", out_tab)

  # Plots
  if (nrow(sum_raw) > 0) {
    sum_raw_plot <- sum_raw %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_sites, decreasing = TRUE)]))

    p_raw <- ggplot(sum_raw_plot, aes(x = SAMPLE, y = ALT_sites)) +
      geom_col(color = "black", fill = "#AA5577") +
      labs(
        title = paste("INDELs (ALT) per sample -", chr),
        x = "Sample",
        y = "Number of ALT INDELs"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_png(p_raw, file.path(plots_dir, paste0("indels_by_sample_", chr, ".png")))
  } else {
    message("No data available for raw count plot in ", chr)
  }

  if (nrow(sum_norm) > 0 && any(!is.na(sum_norm$ALT_per_Mb))) {
    sum_norm_plot <- sum_norm %>%
      filter(!is.na(ALT_per_Mb)) %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_per_Mb, decreasing = TRUE)]))

    if (nrow(sum_norm_plot) > 0) {
      p_norm <- ggplot(sum_norm_plot, aes(x = SAMPLE, y = ALT_per_Mb)) +
        geom_col(color = "black", fill = "#44AA66") +
        labs(
          title = paste("INDELs (ALT) per Mb per sample -", chr),
          x = "Sample",
          y = "INDELs per Mb"
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      save_png(p_norm, file.path(plots_dir, paste0("indels_norm_by_sample_", chr, ".png")))
    } else {
      message("No samples with valid sequence length for normalized plot in ", chr)
    }
  } else {
    message("Insufficient data for normalized plot in ", chr)
  }

  global_counts <- bind_rows(global_counts, sum_raw)
  global_bp     <- bind_rows(global_bp, len_df)
}

# Genome-wide summary per sample
if (nrow(global_counts) > 0) {
  global_raw <- global_counts %>%
    group_by(SAMPLE) %>%
    summarise(ALT_sites = sum(ALT_sites, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ALT_sites))

  global_bp_sum <- global_bp %>%
    group_by(SAMPLE) %>%
    summarise(bp = sum(bp, na.rm = TRUE), .groups = "drop")

  global_norm <- global_raw %>%
    left_join(global_bp_sum, by = "SAMPLE") %>%
    mutate(ALT_per_Mb = if_else(!is.na(bp) & bp > 0, ALT_sites / (bp / 1e6), NA_real_))

  out_global <- file.path(sums_dir, "indels_summary_all_chromosomes_per_sample.csv")
  write_csv(global_norm, out_global)
  message("Genome-wide table saved: ", out_global)

  # Genome-wide raw plot
  if (nrow(global_raw) > 0) {
    global_raw_plot <- global_raw %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_sites, decreasing = TRUE)]))

    p_g1 <- ggplot(global_raw_plot, aes(x = SAMPLE, y = ALT_sites)) +
      geom_col(color = "black", fill = "#AA5577") +
      labs(
        title = "INDELs (ALT) per sample - All chromosomes",
        x = "Sample",
        y = "Number of ALT INDELs"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_png(p_g1, file.path(plots_dir, "indels_by_sample_all.png"))
  }

  # Genome-wide normalized plot
  norm_ok <- global_norm %>% filter(!is.na(ALT_per_Mb))
  if (nrow(norm_ok) > 0) {
    norm_ok <- norm_ok %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_per_Mb, decreasing = TRUE)]))

    p_g2 <- ggplot(norm_ok, aes(x = SAMPLE, y = ALT_per_Mb)) +
      geom_col(color = "black", fill = "#44AA66") +
      labs(
        title = "INDELs per Mb per sample - All chromosomes",
        x = "Sample",
        y = "INDELs per Mb"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_png(p_g2, file.path(plots_dir, "indels_norm_by_sample_all.png"))
  } else {
    message("No samples with valid sequence length for genome-wide normalized plot.")
  }
} else {
  message("No data available for genome-wide summary.")
}
