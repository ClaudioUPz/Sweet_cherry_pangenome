#!/usr/bin/env bash

# =========================================================
# Script: annotate_sv_hits_with_genes.sh
# Description:
# Intersect significant GWAS-associated structural variants
# with gene annotations and generate candidate gene tables.
#
# Main outputs:
# - SVs overlapping genes
# - SVs within 10 kb of genes
# - combined SV-gene annotation table
# - candidate gene ID and gene name lists
# =========================================================

set -euo pipefail

BASE="/path/to/project/8.GWAS"
WORKDIR="${BASE}/mvp_out_mlm_Fru_Crack_St_min20"
GFF="${BASE}/reference_annotation.gff"

cd "$WORKDIR"

#=========================================================
# Step 1: Generate BED file of significant SVs
#=========================================================

awk -F',' 'NR>1 {
    chrom=$2;
    pos=$3;
    sv=$1;
    gsub(/"/,"",chrom);
    gsub(/"/,"",sv);
    print chrom"\t"pos-1"\t"pos"\t"sv
}' Trait.MLM_signals.csv > significant_sv.bed

awk 'BEGIN{OFS="\t"}{
    $1 = "chr"$1;
    print
}' significant_sv.bed > significant_sv.chr.bed

#=========================================================
# Step 2: Generate BED file of genes
#=========================================================

cd "$BASE"

grep -P '\tgene\t' "$GFF" | \
awk -F'\t' 'BEGIN{OFS="\t"}
{
    split($9,a,";");
    id=a[1];
    gsub("ID=","",id);
    print $1,$4,$5,id
}' > genes.bed

#=========================================================
# Step 3: Intersect significant SVs with genes
#=========================================================

cd "$WORKDIR"

bedtools intersect \
  -a significant_sv.chr.bed \
  -b ../genes.bed \
  -wa -wb \
  > significant_sv_in_genes.txt

bedtools window \
  -a significant_sv.chr.bed \
  -b ../genes.bed \
  -w 10000 \
  > significant_sv_near_genes_10kb.txt

#=========================================================
# Step 4: Extract gene ID-to-name table from GFF
#=========================================================

cd "$BASE"

grep -P '\tgene\t' "$GFF" | \
awk -F'\t' 'BEGIN{OFS="\t"}
{
    id="NA";
    name="NA";

    if (match($9,/ID=[^;]+/)) {
        id=substr($9,RSTART+3,RLENGTH-3);
    }
    if (match($9,/Name=[^;]+/)) {
        name=substr($9,RSTART+5,RLENGTH-5);
    }

    print id,name
}' > geneID_to_name.tsv

#=========================================================
# Step 5: Build SV-gene relationship tables
#=========================================================

cd "$WORKDIR"

awk 'BEGIN{OFS="\t"}{print $4,$8,"inside"}' \
  significant_sv_in_genes.txt | sort -u > sv_gene_inside.tsv

awk 'BEGIN{OFS="\t"}{print $4,$8,"near_10kb"}' \
  significant_sv_near_genes_10kb.txt | sort -u > sv_gene_near_10kb.tsv

cat sv_gene_inside.tsv sv_gene_near_10kb.tsv | sort -u > sv_gene_all.tsv

#=========================================================
# Step 6: Add functional gene names
# Final columns:
# SV_ID   gene_ID   gene_name   relation
#=========================================================

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {name[$1]=$2; next}
{
    gene_name = ($2 in name ? name[$2] : "NA");
    print $1,$2,gene_name,$3
}' "$BASE/geneID_to_name.tsv" sv_gene_all.tsv > sv_gene_annotation_all.tsv

#=========================================================
# Step 7: Optional separate outputs
#=========================================================

awk '$4=="inside"' sv_gene_annotation_all.tsv > sv_gene_annotation_inside.tsv
awk '$4=="near_10kb"' sv_gene_annotation_all.tsv > sv_gene_annotation_near_10kb.tsv

#=========================================================
# Step 8: Candidate gene lists
#=========================================================

cut -f2 sv_gene_annotation_all.tsv | sort -u > candidate_geneIDs_10kb.txt
cut -f3 sv_gene_annotation_all.tsv | sort -u > candidate_geneNames_10kb.txt

#=========================================================
# Step 9: Summary
#=========================================================

echo "Annotation completed successfully."
echo "Main output:"
echo "  $WORKDIR/sv_gene_annotation_all.tsv"
echo
echo "Additional outputs:"
echo "  $WORKDIR/sv_gene_annotation_inside.tsv"
echo "  $WORKDIR/sv_gene_annotation_near_10kb.tsv"
echo "  $WORKDIR/candidate_geneIDs_10kb.txt"
echo "  $WORKDIR/candidate_geneNames_10kb.txt"
