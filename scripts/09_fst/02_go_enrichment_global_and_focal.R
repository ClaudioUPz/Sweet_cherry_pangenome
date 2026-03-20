#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Usage:\n",
    "Rscript 02_go_enrichment_global_and_focal.R ",
    "<gene2go.gene.tsv> <genes_in_fst_outliers.list.txt> ",
    "<genes_in_fst_x_gwas.list.txt> <all_genes.list.txt> <output_prefix>\n"
  )
}

gene2go_file <- args[1]
global_file  <- args[2]
focal_file   <- args[3]
background_file <- args[4]
output_prefix <- args[5]

cat("Loading input files...\n")

gene2go <- read_tsv(
  gene2go_file,
  col_names = c("gene", "gos"),
  show_col_types = FALSE
)

required_cols <- c("gene", "gos")
missing_cols <- setdiff(required_cols, colnames(gene2go))
if (length(missing_cols) > 0) {
  stop("Missing required columns in gene2go file: ", paste(missing_cols, collapse = ", "))
}

gene2go_long <- gene2go %>%
  separate_rows(gos, sep = ";") %>%
  filter(str_detect(gos, "^GO:"))

term2gene <- gene2go_long %>%
  select(gos, gene) %>%
  distinct()

genes_with_go <- unique(gene2go_long$gene)

global_genes <- read_lines(global_file)
focal_genes <- read_lines(focal_file)
background_genes <- read_lines(background_file)

global_genes <- intersect(global_genes, genes_with_go)
focal_genes <- intersect(focal_genes, genes_with_go)
background_genes <- intersect(background_genes, genes_with_go)

cat("Genes with GO annotations:", length(genes_with_go), "\n")
cat("Genes in FST outlier regions:", length(global_genes), "\n")
cat("Genes in FST-GWAS overlapping regions:", length(focal_genes), "\n")
cat("Background genes with GO annotations:", length(background_genes), "\n\n")

run_enrichment <- function(gene_list, label, universe_genes, term2gene_table, output_prefix) {
  if (length(gene_list) < 5) {
    cat(label, "skipped: too few genes\n")
    return(NULL)
  }

  ego <- enricher(
    gene = gene_list,
    universe = universe_genes,
    TERM2GENE = term2gene_table,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  if (is.null(ego)) {
    cat(label, "no enriched terms detected\n")
    return(NULL)
  }

  res <- as.data.frame(ego)

  if (nrow(res) == 0) {
    cat(label, "no significant GO terms\n")
    return(NULL)
  }

  output_tsv <- paste0(output_prefix, ".", label, ".tsv")
  write_tsv(res, output_tsv)
  cat("Saved:", output_tsv, "\n")

  png(
    paste0(output_prefix, ".", label, ".dotplot.png"),
    width = 1600,
    height = 1000,
    res = 160
  )
  print(dotplot(ego, showCategory = 20))
  dev.off()

  pdf(
    paste0(output_prefix, ".", label, ".dotplot.pdf"),
    width = 8,
    height = 5
  )
  print(dotplot(ego, showCategory = 20))
  dev.off()

  png(
    paste0(output_prefix, ".", label, ".barplot.png"),
    width = 1600,
    height = 1000,
    res = 160
  )
  print(barplot(ego, showCategory = 20))
  dev.off()

  return(res)
}

cat("Running enrichment for genes in FST outlier regions...\n")
global_results <- run_enrichment(
  gene_list = global_genes,
  label = "global",
  universe_genes = background_genes,
  term2gene_table = term2gene,
  output_prefix = output_prefix
)

cat("\nRunning enrichment for genes in FST-GWAS overlapping regions...\n")
focal_results <- run_enrichment(
  gene_list = focal_genes,
  label = "focal",
  universe_genes = background_genes,
  term2gene_table = term2gene,
  output_prefix = output_prefix
)
