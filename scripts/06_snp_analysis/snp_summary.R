# 06_summary_per_sample_from_step2.R
# Resumen por MUESTRA (sin grupos) a partir de los CSV del Paso 2:
# - Conteo bruto de SNPs ALT por muestra y cromosoma
# - Normalización por longitud de cada muestra (bp) en ese cromosoma, leída de 5.pangenome/chrN/<SAMPLE>_chrN.fa
# - Plots robustos (device png) y tablas por cromosoma y globales

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr)
  library(ggplot2); library(viridis); library(tidyr); library(stringr)
})
options(dplyr.summarise.inform = FALSE)

# ------------ CONFIG ------------
BASE <- "/path/to/basedir"
chr_ids <- paste0("chr", 1:8)      # para probar solo chr1: chr_ids <- "chr1"
details_name <- "snps_details_%s.csv"    # "%s" -> chr1, chr2, ...
plots_dir <- file.path(BASE, "SNP_dir", "plots_per_sample")
sums_dir  <- file.path(BASE, "SNP_dir", "summaries_per_sample")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sums_dir,  showWarnings = FALSE, recursive = TRUE)

# PNG robusto (evita PNG en blanco en headless)
save_png <- function(plot, file, width=12, height=8, dpi=300, bg="white"){
  png(filename = file, width = width*dpi, height = height*dpi, res = dpi, bg = bg)
  print(plot); dev.off()
  message("PNG guardado: ", file)
}

# Longitud total (bp) de un FASTA (ignora líneas ">" de cabecera)
fasta_bp <- function(path){
  if (!file.exists(path)) return(NA_real_)
  L <- readLines(path, warn = FALSE)
  sum(nchar(L[!startsWith(L, ">")]))
}

# ------------ ACUMULADORES ------------
global_counts <- tibble(SAMPLE=character(), ALT_sites=double())
global_bp     <- tibble(SAMPLE=character(), bp=double())

# ------------ LOOP POR CROMOSOMA ------------
for (chr in chr_ids) {
  message("\n========== ", chr, " ==========")
  csv <- file.path(BASE, "SNP_dir", chr, sprintf(details_name, chr))
  smp <- file.path(BASE, "SNP_dir", chr, "samples.txt")

  if (!file.exists(csv)) { message("Saltando (no existe): ", csv); next }

  # CSV del Paso 2: CHROM,POS,REF,ALT,SAMPLE,GT
  d <- read_csv(csv, show_col_types = FALSE,
                col_types = cols(
                  CHROM  = col_character(),
                  POS    = col_double(),
                  REF    = col_character(),
                  ALT    = col_character(),
                  SAMPLE = col_character(),
                  GT     = col_character()
                ))
  message("Filas ALT en ", chr, ": ", format(nrow(d), big.mark=","))

  # Muestras presentes en el VCF de este cromosoma
  samples_chr <- if (file.exists(smp)) readLines(smp) else sort(unique(d$SAMPLE))

  # Conteo bruto por muestra
  sum_raw <- d %>%
    count(SAMPLE, name = "ALT_sites") %>%
    arrange(desc(ALT_sites))

  # Longitud por muestra en este cromosoma (FASTA: 5.pangenome/chrN/<SAMPLE>_chrN.fa)
  fasta_paths <- file.path(BASE, "pangenome_dir", chr, paste0(samples_chr, "_", chr, ".fa"))
  bp_vec <- map_dbl(fasta_paths, fasta_bp)
  len_df <- tibble(SAMPLE = samples_chr, bp = bp_vec)

  # Join y normalización por Mb de cada muestra
  sum_norm <- sum_raw %>%
    left_join(len_df, by = "SAMPLE") %>%
    mutate(ALT_per_Mb = if_else(!is.na(bp) & bp > 0, ALT_sites / (bp/1e6), NA_real_))

  # Guardar tabla por cromosoma
  out_tab <- file.path(sums_dir, paste0("snps_summary_", chr, "_per_sample.csv"))
  write_csv(sum_norm, out_tab)
  message("Tabla guardada: ", out_tab)

  # ---- PLOTS ----
  if (nrow(sum_raw) > 0) {
    sum_raw_plot <- sum_raw %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_sites, decreasing = TRUE)]))
    p_raw <- ggplot(sum_raw_plot, aes(x = SAMPLE, y = ALT_sites)) +
      geom_col(color = "black", fill = "#5B002C") +  # burdeo oscuro
      labs(title = paste("SNPs (ALT) per sample -", chr),
           x = "Sample", y = "N° SNPs (ALT)") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "#F6ECEF", color = NA),
        plot.background  = element_rect(fill = "#F6ECEF", color = NA),
        text = element_text(color = "#333333"),
        axis.title = element_text(face = "bold", color = "#5B002C"),
        plot.title = element_text(face = "bold", color = "#5B002C")
      )
    save_png(p_raw, file.path(plots_dir, paste0("snps_by_sample_burdeo_", chr, ".png")))
    } else message("Sin datos para gráfico bruto en ", chr)

  if (nrow(sum_norm) > 0 && any(!is.na(sum_norm$ALT_per_Mb))) {
    sum_norm_plot <- sum_norm %>%
      filter(!is.na(ALT_per_Mb)) %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_per_Mb, decreasing = TRUE)]))
    if (nrow(sum_norm_plot) > 0) {
      p_norm <- ggplot(sum_norm_plot, aes(x = SAMPLE, y = ALT_per_Mb)) +
        geom_col(color = "black", fill = "#A8334C") +  # granate medio (Burdeo Regina)
        labs(title = paste("SNPs/Mb -", chr),
             x = "Sample", y = "SNPs ALT per Mb") +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "#F6ECEF", color = NA),
          plot.background  = element_rect(fill = "#F6ECEF", color = NA),
          text = element_text(color = "#333333"),
          axis.title = element_text(face = "bold", color = "#5B002C"),
          plot.title = element_text(face = "bold", color = "#5B002C")
        )
      save_png(p_norm, file.path(plots_dir, paste0("snps_norm_by_sample_burdeo_", chr, ".png")))
    } else {
      message("No hay muestras con longitud válida para el gráfico normalizado en ", chr)
    }
  } else message("Sin datos suficientes para gráfico normalizado en ", chr)
  
  # Acumular para global
  global_counts <- bind_rows(global_counts, sum_raw)
  global_bp     <- bind_rows(global_bp, len_df)
}

# ------------ RESUMEN GLOBAL POR MUESTRA ------------
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
    mutate(ALT_per_Mb = if_else(!is.na(bp) & bp > 0, ALT_sites / (bp/1e6), NA_real_))
  
  out_global <- file.path(sums_dir, "snps_summary_all_chromosomes_per_sample.csv")
  write_csv(global_norm, out_global)
  message("Tabla global guardada: ", out_global)
  
  # Plot global bruto
  if (nrow(global_raw) > 0) {
    global_raw_plot <- global_raw %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_sites, decreasing = TRUE)]))
    p_g1 <- ggplot(global_raw_plot, aes(x = SAMPLE, y = ALT_sites)) +
      geom_col(color = "black", fill = "#5B002C") +  # burdeo oscuro
      labs(title = "SNPs (ALT) per sample - All chromosomes",
           x = "Muestra", y = "N° SNPs (ALT)") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "#F6ECEF", color = NA),
        plot.background  = element_rect(fill = "#F6ECEF", color = NA),
        text = element_text(color = "#333333"),
        axis.title = element_text(face = "bold", color = "#5B002C"),
        plot.title = element_text(face = "bold", color = "#5B002C")
      )
    save_png(p_g1, file.path(plots_dir, "snps_by_sample_all_burdeo.png"))
  }
  
  # Plot global normalizado
  norm_ok <- global_norm %>% filter(!is.na(ALT_per_Mb))
  if (nrow(norm_ok) > 0) {
    norm_ok <- norm_ok %>%
      mutate(SAMPLE = factor(SAMPLE, levels = SAMPLE[order(ALT_per_Mb, decreasing = TRUE)]))
    p_g2 <- ggplot(norm_ok, aes(x = SAMPLE, y = ALT_per_Mb)) +
      geom_col(color = "black", fill = "#A8334C") +  # granate medio
      labs(title = "SNPs/Mb - All chromosomes",
           x = "Sample", y = "SNPs/Mb") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "#F6ECEF", color = NA),
        plot.background  = element_rect(fill = "#F6ECEF", color = NA),
        text = element_text(color = "#333333"),
        axis.title = element_text(face = "bold", color = "#5B002C"),
        plot.title = element_text(face = "bold", color = "#5B002C")
      )
    save_png(p_g2, file.path(plots_dir, "snps_norm_by_sample_all_burdeo.png"))
  } else {
    message("No hay muestras con longitud válida para el gráfico global normalizado.")
  }
  
} else {
  message("No hubo datos para resumen global (¿faltan CSV del Paso 2?).")
}
