# 09_Fst – Population differentiation and functional analysis

This directory contains the workflow to compute population differentiation (FST) and integrate it with GWAS signals, gene annotation, and functional enrichment analyses.

The pipeline identifies genomic regions under selection and links them to candidate genes and biological processes.

---

## Overview of the workflow

The FST analysis is divided into three main components:

### 1. FST computation
- SNP-level FST calculation (Hudson estimator)
- Aggregation into genomic windows
- Identification of outlier regions

### 2. Genomic interpretation
- Extraction of genes in FST outlier regions
- Integration with GWAS signals
- Identification of candidate genes

### 3. Functional analysis
- GO enrichment analysis
- Mapping GO terms to genomic regions
- Visualization of chromosome contributions

---

## Requirements

### Software
- bcftools
- samtools
- bedtools
- awk, sort, cut, grep

### R packages
- dplyr
- readr
- ggplot2
- tidyr
- stringr
- clusterProfiler
- GO.db

---

## Input data

- SNP VCF (filtered cohort)
- Sample group file (sample_groups_checked.tsv)
- Gene annotation (.gff3)
- GWAS results (Trait.GLM_signals.csv)
- InterPro annotation (*_vs_interpro.tsv)

---

## Workflow

---

### STEP 13 – Plot FST values for all chromosomes

```bash
conda activate r_fst_env

mkdir -p /path/to/project/selection_vcf/Fst_plots
cd /path/to/project/selection_vcf/Fst_plots

Rscript /path/to/project/scripts/plot_fst_with_outliers_all_chr.R \
  /path/to/project/selection_vcf/snps.Fst_landrace_vs_modern.windows_100kb.tsv \
  /path/to/project/selection_vcf/Fst_outliers \
  Fst_outliers

Output:

One plot per chromosome highlighting FST outlier windows

STEP 14 – Compare chromosome-level FST patterns
Rscript /path/to/project/scripts/compare_chr_selection_fst.R \
  /path/to/project/selection_vcf/snps.Fst_landrace_vs_modern.windows_100kb.tsv \
  /path/to/project/selection_vcf/Fst_genomewide_p99

Output:

Chromosome summary table

Barplots:

number of top windows

proportion of top windows

maximum FST per chromosome

Downstream analysis: genes and functional annotation
STEP 15 – Extract genes from FST regions and build gene2GO
bash /path/to/project/scripts/01_fst_regions_genes_and_gene2go.sh \
  /path/to/project/selection_vcf/Fst_outliers \
  /path/to/project/reference/reference_annotation.gff3 \
  /path/to/project/selection_vcf/selection_genes \
  100000 \
  /path/to/project/gwas/Trait.GLM_signals.csv \
  250000 \
  /path/to/project/annotation/reference_vs_interpro.tsv \
  /path/to/project/gwas/gwas_sv_hits.bed

Output:

genes_in_fst_outliers.list.txt

genes_in_fst_x_gwas.list.txt

all_genes.list.txt

gene2go.gene.tsv

merged FST regions (BED)

STEP 16 – GO enrichment analysis
Rscript /path/to/project/scripts/02_go_enrichment_global_and_focal.R \
  /path/to/project/selection_vcf/selection_genes/gene2go.gene.tsv \
  /path/to/project/selection_vcf/selection_genes/genes_in_fst_outliers.list.txt \
  /path/to/project/selection_vcf/selection_genes/genes_in_fst_x_gwas.list.txt \
  /path/to/project/selection_vcf/selection_genes/all_genes.list.txt \
  /path/to/project/selection_vcf/selection_genes/GO

Output:

GO enrichment tables (global and focal)

Dotplots and barplots

STEP 17 – Add GO term names
Rscript /path/to/project/scripts/03_add_go_names.R \
  /path/to/project/selection_vcf/selection_genes/GO.global.tsv \
  /path/to/project/selection_vcf/selection_genes/GO.global.named.tsv
STEP 18 – Trace genomic origin of GO terms
bash /path/to/project/scripts/05_trace_go_origin.sh \
  /path/to/project/selection_vcf/selection_genes/GO.global.tsv \
  /path/to/project/selection_vcf/selection_genes/genes.bed \
  /path/to/project/selection_vcf/selection_genes/fst_outliers_merged_d100000.chromfix.bed \
  /path/to/project/selection_vcf/selection_genes/GO_trace

Output:

GO-to-gene mapping

Chromosome contribution per GO

Genomic span of GO-associated genes

Summary table linking GO terms to FST regions

STEP 19 – Plot chromosome contribution to GO terms
Rscript /path/to/project/scripts/06_plot_go_chr_dominance.R \
  /path/to/project/selection_vcf/selection_genes/GO_trace/go_by_chr.tsv \
  /path/to/project/selection_vcf/selection_genes/GO.global.named.tsv \
  /path/to/project/selection_vcf/selection_genes/GO_trace/go_chr_dominance_named

Output:

Stacked barplots showing chromosome contribution per GO term

STEP 20 – Extract windows contributing to GO terms
Rscript /path/to/project/scripts/07_go_windows_positions.R \
  /path/to/project/selection_vcf/selection_genes/GO_trace/go_genes_in_outlier_windows.tsv \
  /path/to/project/selection_vcf/selection_genes/GO_trace/go_windows

Output:

Table of FST windows contributing to each GO term

Visualization of genomic positions

STEP 21 – Generate GO gene table (supplementary table)
bash /path/to/project/scripts/08_make_go_gene_table.sh \
  /path/to/project/selection_vcf/selection_genes/GO.global.tsv \
  /path/to/project/selection_vcf/selection_genes/GO_trace/go_to_gene.coords.tsv \
  /path/to/project/selection_vcf/selection_genes/GO_trace/supplementary_table_go_genes.tsv

Output:

Table linking GO terms, genes, and genomic coordinates

Key concepts

FST (Hudson estimator): measures genetic differentiation between populations

Outliers: defined as top FST windows (per chromosome or genome-wide threshold)

FST–GWAS integration: identifies regions under selection associated with phenotypic traits

GO enrichment: identifies biological processes overrepresented in selected regions

Notes

Chromosome names are normalized (e.g., chr4__Regina → chr4)

Window size is configurable (default: 100 kb)

GWAS integration includes SNP and optional SV signals

GO enrichment is based on custom gene2go mapping derived from InterPro
