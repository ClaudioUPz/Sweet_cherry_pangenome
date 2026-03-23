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
- plink
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

## Linkage Disequilibrium (LD) decay analysis

To support the selection of the window size used for FST analyses, genome-wide linkage disequilibrium (LD) decay was estimated using PLINK v1.9.

### Step 1: LD computation

```bash
bash run_plink_LD.sh
Step 2: LD decay summarization
nohup ./run_ld_decay_parallel.sh > logs/ld_decay_parallel.log 2>&1 &

Output:

output/plink_ld/LD_decay_summary.tsv
Step 3: LD decay visualization
Rscript plot_LD_decay.R \
  output/plink_ld/LD_decay_summary.tsv \
  output/plink_ld/LD_decay

Generated outputs:

LD_decay.png / LD_decay.pdf
LD_decay_zoom200kb.png / LD_decay_zoom200kb.pdf
LD_decay_summary.txt
Interpretation

LD decayed rapidly within the first ~100 kb, reaching r² ≈ 0.2 at ~100 kb, and decreased more gradually at larger distances.

The chosen window size (100 kb) corresponds to the distance at which LD decays to r² ≈ 0.2, supporting its use in FST analyses.

All LD analyses were performed using the same filtered SNP dataset used for FST estimation (MAF ≥ 0.01, missingness < 0.2).

Workflow
STEP 13 – Plot FST values for all chromosomes
conda activate r_fst_env

mkdir -p output/Fst_plots
cd output/Fst_plots

Rscript scripts/plot_fst_with_outliers_all_chr.R \
  input/snps.Fst_landrace_vs_modern.windows_100kb.tsv \
  input/Fst_outliers \
  Fst_outliers

Output:

One plot per chromosome highlighting FST outlier windows
STEP 14 – Compare chromosome-level FST patterns
Rscript scripts/compare_chr_selection_fst.R \
  input/snps.Fst_landrace_vs_modern.windows_100kb.tsv \
  input/Fst_genomewide_p99

Output:

Chromosome summary table
Barplots:
number of top windows
proportion of top windows
maximum FST per chromosome
STEP 15 – Extract genes from FST regions and build gene2GO
bash scripts/01_fst_regions_genes_and_gene2go.sh \
  input/Fst_outliers \
  input/reference_annotation.gff3 \
  output/selection_genes \
  100000 \
  input/Trait.GLM_signals.csv \
  250000 \
  input/reference_vs_interpro.tsv \
  input/gwas_sv_hits.bed

Output:

genes_in_fst_outliers.list.txt
genes_in_fst_x_gwas.list.txt
all_genes.list.txt
gene2go.gene.tsv
merged FST regions (BED)
STEP 16 – GO enrichment analysis
Rscript scripts/02_go_enrichment_global_and_focal.R \
  output/selection_genes/gene2go.gene.tsv \
  output/selection_genes/genes_in_fst_outliers.list.txt \
  output/selection_genes/genes_in_fst_x_gwas.list.txt \
  output/selection_genes/all_genes.list.txt \
  output/selection_genes/GO

Output:

GO enrichment tables
Dotplots and barplots
STEP 17 – Add GO term names
Rscript scripts/03_add_go_names.R \
  output/selection_genes/GO.global.tsv \
  output/selection_genes/GO.global.named.tsv
STEP 18 – Trace genomic origin of GO terms
bash scripts/05_trace_go_origin.sh \
  output/selection_genes/GO.global.tsv \
  output/selection_genes/genes.bed \
  output/selection_genes/fst_outliers_merged_d100000.chromfix.bed \
  output/selection_genes/GO_trace

Output:

GO-to-gene mapping
Chromosome contribution per GO
Genomic span of GO-associated genes
STEP 19 – Plot chromosome contribution to GO terms
Rscript scripts/06_plot_go_chr_dominance.R \
  output/selection_genes/GO_trace/go_by_chr.tsv \
  output/selection_genes/GO.global.named.tsv \
  output/selection_genes/GO_trace/go_chr_dominance_named
STEP 20 – Extract windows contributing to GO terms
Rscript scripts/07_go_windows_positions.R \
  output/selection_genes/GO_trace/go_genes_in_outlier_windows.tsv \
  output/selection_genes/GO_trace/go_windows
STEP 21 – Generate GO gene table (supplementary table)
bash scripts/08_make_go_gene_table.sh \
  output/selection_genes/GO.global.tsv \
  output/selection_genes/GO_trace/go_to_gene.coords.tsv \
  output/selection_genes/GO_trace/supplementary_table_go_genes.tsv
Key concepts
FST (Hudson estimator): measures genetic differentiation between populations
Outliers: top FST windows
FST–GWAS integration: links selection signals to traits
GO enrichment: identifies biological processes under selection
LD decay: informs appropriate genomic window size
Notes
Chromosome names are normalized (e.g., chr4__Regina → chr4)
Window size is configurable (default: 100 kb)
GWAS integration includes SNP and optional SV signals
GO enrichment is based on gene2GO derived from InterPro
The selected FST window size (100 kb) is supported by LD decay analysis
