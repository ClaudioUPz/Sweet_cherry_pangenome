# 🍒 Sweet Cherry Pangenome Analysis

This repository contains the complete computational workflow for the construction, analysis, and biological interpretation of a **sweet cherry (Prunus avium) pangenome**, including:

- Genome assembly and integration
- Structural variant (SV) detection
- SNP-based analyses
- Genome-wide association studies (GWAS)
- Population differentiation (FST)
- Functional annotation and GO enrichment

The pipeline integrates **pangenome variation, GWAS, and selection signals** to identify candidate genomic regions underlying agronomic traits.

---

## 📁 Repository structure

```bash
Sweet_cherry_pangenome/
│
├── 01_preprocessing/
├── 02_pangenome/
├── 03_alignment/
├── 04_variant_calling/
├── 05_structural_variants/
├── 06_annotation/
├── 07_population_genomics/
├── 08_GWAS/
├── 09_Fst/
│
└── README.md
🔬 Conceptual workflow
Raw reads
   ↓
Genome assemblies
   ↓
Pangenome graph
   ↓
Variant detection (SNPs + SVs)
   ↓
├── GWAS (trait associations)
└── FST (selection signals)
        ↓
Integration (GWAS × FST × annotation)
        ↓
Candidate genes
        ↓
GO enrichment & biological interpretation
⚙️ Requirements
Software

bcftools

samtools

bedtools

minimap2 / vg (upstream steps)

awk, sort, grep, cut

R packages

dplyr

readr

ggplot2

bigmemory

rMVP

tidyr

stringr

clusterProfiler

GO.db

📊 Modules
08_GWAS – Genome-wide association studies

This module performs GWAS using both:

Structural variants (SVs)

SNPs

Key features

Joint SNP + SV analysis

rMVP (GLM and MLM models)

Population structure correction (PCA + kinship)

Manhattan plots and signal extraction

Outputs

Significant loci

Manhattan plots

Candidate variant lists

Regions for downstream annotation

09_Fst – Population differentiation and selection

This module identifies genomic regions under selection and links them to biological functions.

🔹 FST computation

SNP-level FST (Hudson estimator)

Aggregation into genomic windows (e.g., 100 kb)

Identification of outlier regions

🔹 Integration with genomic features

Intersection with gene annotation

Integration with GWAS signals

Identification of candidate genes under selection

🔹 Functional analysis

GO enrichment (clusterProfiler)

Mapping GO terms to genomic regions

Chromosome-level contribution analysis

🧪 FST workflow summary

Sample group definition (landrace, early selection, modern)

VCF filtering (SNP + SV)

Genotype matrix extraction

SNP-level FST calculation

Window-based aggregation

Outlier detection

Visualization

Gene extraction

GWAS integration

GO enrichment

Chromosome contribution analysis

📈 Key outputs
FST

SNP-level FST values

Window-based FST tables

Outlier regions (top 5% or genome-wide thresholds)

Integration

Genes in FST regions

Genes overlapping GWAS signals

Candidate gene lists

Functional

GO enrichment tables

GO–gene mappings

Chromosome contribution plots

Supplementary tables for publication

🔗 Integration across analyses

This repository is designed to connect multiple layers of genomic information:

Layer	Description
Pangenome	Captures structural diversity
SNPs	Fine-scale variation
SVs	Large-effect variants
GWAS	Trait associations
FST	Selection signals
GO	Functional interpretation
🧠 Key concepts

Pangenome: graph-based representation of multiple genomes

GWAS: association between variants and traits

FST (Hudson): genetic differentiation between populations

Outliers: genomic regions under potential selection

GO enrichment: overrepresented biological processes

📌 Notes

Chromosome naming is normalized (e.g., chr4__Regina → chr4)

Window size is configurable (default: 100 kb)

GWAS integrates SNPs and structural variants

FST is computed between landrace vs modern breeding

Functional annotation is based on InterPro-derived mappings

🚀 Future extensions

Integration with 3D genome (Hi-C / Micro-C)

Structural variant functional impact analysis

Pangenome-aware GWAS

Cross-species comparative genomics

## 👤 Maintainer

Claudio Urra Pérez  
Postdoctoral researcher – Plant Genomics  

---

## 📄 Project context

This repository is part of a collaborative research project on the sweet cherry pangenome.

For full authorship, contributions, and project leadership, please refer to the associated manuscript.
📄 Citation

If you use this pipeline, please cite the corresponding pangenome and GWAS manuscripts (in preparation / submitted).
