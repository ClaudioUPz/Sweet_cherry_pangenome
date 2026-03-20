#!/usr/bin/env bash

set -euo pipefail

module load bcftools 2>/dev/null || true

if [[ $# -ne 3 ]]; then
  echo "Usage:"
  echo "  bash export_snp_genotypes_for_fst.sh <filtered_snp_vcf.gz> <genotypes_output.tsv> <sample_order_output.txt>"
  exit 1
fi

INPUT_VCF=$1
GENO_OUT=$2
SAMPLE_ORDER_OUT=$3

echo "Input SNP VCF:        $INPUT_VCF"
echo "Genotype table out:   $GENO_OUT"
echo "Sample order file:    $SAMPLE_ORDER_OUT"
echo

bcftools query \
  -H \
  -f '%CHROM\t%POS[\t%GT]\n' \
  "$INPUT_VCF" \
  > "$GENO_OUT"

bcftools query -l "$INPUT_VCF" > "$SAMPLE_ORDER_OUT"

echo "Genotype table written to: $GENO_OUT"
echo "Sample order file written to: $SAMPLE_ORDER_OUT"
