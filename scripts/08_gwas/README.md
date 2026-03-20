# Step 8: GWAS analysis

This step describes the preparation of short-read data, read mapping, variant calling, and genome-wide association analyses performed using both structural variants (SVs) and SNPs.

---

## Overview

The workflow includes:

1. Preparing sample lists for parallel processing
2. Mapping short reads to graph-based and linear references
3. Calling structural variants from graph-based alignments
4. Calling SNPs and exporting genotype matrices
5. Running GWAS with structural variants
6. Running GWAS with SNPs

---

## Requirements

- bash
- SLURM
- samtools
- vg
- minimap2
- bcftools
- delly
- R
- Conda
- Required R packages:
  - bigmemory
  - rMVP
  - readr
  - dplyr
  - readxl
  - ggplot2
  - ragg

---

## Step 1: Prepare sample lists

Before mapping, ordered sample lists were generated to allow parallel processing of paired-end short-read datasets.

### Script

- `prepare_sample_lists.sh`

### Input

- Sample directories under:

```text
/path/to/project/2.short_reads/

Each sample directory is expected to contain paired-end FASTQ files with _1 and _2 in the filename.

Output

ID_list.txt

R1_list.txt

R2_list.txt

These files contain synchronized sample IDs and paths to forward and reverse reads.

Step 2: Read mapping

Short reads were mapped against:

a graph-based pangenome using VG Giraffe

a linear reference using Minimap2

This enabled comparison between graph-based and linear mapping strategies.

Script

map_reads_graph_vs_linear.sh

Input

ID_list_part_X.txt

R1_list_part_X.txt

R2_list_part_X.txt

merged graph in GFA/GBZ format

combined linear reference FASTA

Output

graph-based BAM files

linear BAM files

mapping summary table

Step 3: Structural variant calling

Structural variants were called from graph-based mapped reads.

The workflow included:

conversion of graph alignments to linear coordinates using vg inject and vg surject

BAM sorting and indexing

merging replicate BAM files

DELLY SV calling per sample

construction of union SV sites

joint genotyping across all samples

generation of filtered multi-sample SV VCF files

export of genotype matrices for GWAS

Script

call_sv_delly_from_graph.sh

Main outputs

delly.merged.vcf.gz

delly.merged.split.filled.vcf.gz

filtered SV panels:

delly.SV_min20.noBND.vcf.gz

delly.SV_pass_precise_min50.noBND.vcf.gz

rMVP input files:

sv_matrix_num_forGWAS_min20.txt

sv_map_forGWAS_min20.txt

sv_map_forGWAS_min20.numeric.txt

Step 4: SNP calling and export

SNPs were called from merged BAM files using bcftools mpileup and bcftools call, followed by filtering and export to rMVP-compatible genotype and map files.

Script

call_snps_for_gwas.sh

Filtering criteria

minimum mapping quality: 20

minimum base quality: 20

SNPs only

QUAL >= 30

F_MISSING < 0.2

AF >= 0.01

AF <= 0.99

Main outputs

snp_matrix_num.txt

snp_map.txt

snp_map.numeric.txt

id_mapping_snps.csv

Step 5: GWAS with structural variants

GWAS with structural variants was performed using rMVP under both GLM and MLM models.

The workflow included:

phenotype filtering

genotype/phenotype sample matching

marker filtering by call rate and minor allele frequency

imputation of missing genotypes by marker mean

kinship estimation with VanRaden method

principal component analysis

association testing using GLM and MLM

Manhattan plot generation

export of significant hits and region files

phenotype and QQ diagnostics

Script

run_gwas_sv_rMVP.R

Input files

phenotype file:

Description_accessions_short_read_and_BLUPs.xlsx

SV genotype matrix:

sv_matrix_num_forGWAS_min20.txt

SV map files:

sv_map_forGWAS_min20.txt

sv_map_forGWAS_min20.numeric.txt

Main outputs

GLM results directory:

mvp_out_glm_<TRAIT>_min20/

MLM results directory:

mvp_out_mlm_<TRAIT>_min20/

Manhattan plots:

plots_gwas_png/Manhattan_GLM_<TRAIT>_SV.png

plots_gwas_png/Manhattan_MLM_<TRAIT>_SV.png

annotated hits:

hits_<TRAIT>_min20_annot.csv

bcftools region file:

hits_<TRAIT>_min20.regions.tsv

Step 6: GWAS with SNPs

GWAS with SNPs was performed using rMVP under both GLM and MLM models, following the same general workflow used for structural variants.

Script

run_gwas_snp_rMVP.R

Input files

phenotype file:

Description_accessions_short_read_and_BLUPs.xlsx

SNP genotype matrix:

snp_matrix_num.txt

SNP map files:

snp_map.txt

snp_map.numeric.txt

Main outputs

GLM results directory:

mvp_out_glm_<TRAIT>_SNP_min20/

MLM results directory:

mvp_out_mlm_<TRAIT>_SNP_min20/

Manhattan plots:

plots_gwas_png/Manhattan_GLM_<TRAIT>_SNP.png

plots_gwas_png/Manhattan_MLM_<TRAIT>_SNP.png

annotated hits:

hits_<TRAIT>_SNP_min20_annot.csv

bcftools region file:

hits_<TRAIT>_SNP_min20.regions.tsv

Notes

Phenotypic BLUP values were read directly from the input spreadsheet and were not recalculated in the GWAS scripts.

Sample identifiers from genotype matrices were simplified before matching them with phenotype accession codes.

Marker filtering was applied before association testing:

call rate ≥ 0.8

minor allele frequency ≥ 0.01

Missing genotypes were imputed using the mean genotype value per marker.

Both GLM and MLM models were run for each trait.

The same GWAS framework was used for both SV and SNP datasets, enabling direct comparison of association results.
