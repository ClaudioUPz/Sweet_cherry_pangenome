#!/usr/bin/env bash

set -euo pipefail

module load bcftools 2>/dev/null || true
module load samtools 2>/dev/null || true

if [[ $# -ne 4 ]]; then
  echo "Usage:"
  echo "  bash filter_vcfs_for_selection.sh <sample_groups_checked.tsv> <snps.vcf.gz> <sv.vcf.gz> <output_dir>"
  exit 1
fi

SAMPLE_GROUPS_TSV=$1
VCF_SNP=$2
VCF_SV=$3
OUTDIR=$4

echo "Sample groups table: $SAMPLE_GROUPS_TSV"
echo "SNP VCF:             $VCF_SNP"
echo "SV VCF:              $VCF_SV"
echo "Output directory:    $OUTDIR"
echo

mkdir -p "$OUTDIR"

echo "Generating sample list from group table..."
SAMPLE_LIST="$OUTDIR/samples_for_selection.txt"

cut -f1 "$SAMPLE_GROUPS_TSV" | tail -n +2 > "$SAMPLE_LIST"

N_SAMPLES=$(wc -l < "$SAMPLE_LIST")
echo "Number of samples retained for selection analysis: $N_SAMPLES"
echo "Example sample IDs:"
head "$SAMPLE_LIST"
echo

echo "Filtering SNP VCF..."
SNP_OUT="$OUTDIR/snps.selection_samples.vcf.gz"

bcftools view \
  -S "$SAMPLE_LIST" \
  -Oz -o "$SNP_OUT" \
  "$VCF_SNP"

bcftools index -t "$SNP_OUT"

echo "Filtered SNP VCF written to: $SNP_OUT"
echo

echo "Filtering SV VCF..."
SV_OUT="$OUTDIR/sv.selection_samples.vcf.gz"

bcftools view \
  -S "$SAMPLE_LIST" \
  -Oz -o "$SV_OUT" \
  "$VCF_SV"

bcftools index -t "$SV_OUT"

echo "Filtered SV VCF written to: $SV_OUT"
echo

echo "Selection VCF preparation completed:"
echo " - Sample list:      $SAMPLE_LIST"
echo " - Filtered SNP VCF: $SNP_OUT"
echo " - Filtered SV VCF:  $SV_OUT"
