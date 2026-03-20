#!/usr/bin/env bash

set -euo pipefail

module load bcftools 2>/dev/null || true
module load samtools 2>/dev/null || true

if [[ $# -ne 4 ]]; then
  echo "Usage:"
  echo "  bash prepare_sample_groups.sh <sample_groups_from_blup.tsv> <snps.vcf.gz> <sv.vcf.gz> <output.tsv>"
  exit 1
fi

BLUP_TSV=$1
VCF_SNP=$2
VCF_SV=$3
OUTFILE=$4

echo "BLUP table: $BLUP_TSV"
echo "SNP VCF:    $VCF_SNP"
echo "SV VCF:     $VCF_SV"
echo "Output:     $OUTFILE"
echo

echo "Extracting sample names from SNP VCF..."
bcftools query -l "$VCF_SNP" | sort > snp_samples.txt

echo "Extracting sample names from SV VCF..."
bcftools query -l "$VCF_SV" | sort > sv_samples.txt

echo "Extracting sample names from BLUP group table..."
cut -f1 "$BLUP_TSV" | tail -n +2 | sort > blup_samples.txt

echo "Comparing sample lists..."

comm -12 blup_samples.txt snp_samples.txt > tmp_match_snps.txt
comm -12 blup_samples.txt sv_samples.txt  > tmp_match_svs.txt
comm -12 tmp_match_snps.txt tmp_match_svs.txt > good_samples.txt

echo -e "SampleID\tStatus" > "$OUTFILE"

awk 'NR==FNR {status[$1]=$2; next}
     ($1 in status) {print $1 "\t" status[$1]}' \
    "$BLUP_TSV" \
    good_samples.txt >> "$OUTFILE"

echo
echo "Output written to: $OUTFILE"
echo "Samples shared across BLUP, SNP, and SV datasets: $(wc -l < good_samples.txt)"
echo

echo "Diagnostic summary"
echo "------------------"

echo "BLUP samples missing from SNP VCF:"
comm -23 blup_samples.txt snp_samples.txt || true

echo
echo "BLUP samples missing from SV VCF:"
comm -23 blup_samples.txt sv_samples.txt || true

echo
echo "SNP samples missing from SV VCF:"
comm -23 snp_samples.txt sv_samples.txt || true

echo
echo "Done."
