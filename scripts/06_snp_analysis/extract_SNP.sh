#!/bin/bash
#SBATCH --job-name=snps_vcf
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -o logs/snp_vcf.%j.out
#SBATCH -e logs/snp_vcf.%j.err

set -euo pipefail

# Load module
module load bcftools

BASE="/path/to/basedir"
IN_DIR="$BASE/pangenome_dir"
OUT_DIR="$BASE/SNP_dir"

mkdir -p "$OUT_DIR" logs

for chr in {1..8}; do
  in_vcf="$IN_DIR/chr${chr}/output/chr${chr}_run.vcf.gz"
  if [[ ! -f "$in_vcf" ]]; then
    echo "ERROR: Doesn't exist $in_vcf"; exit 1
  fi

  out_chr_dir="$OUT_DIR/chr${chr}"
  mkdir -p "$out_chr_dir"

  # SNPs (keep bi y multi-allelic). Index with tabix.
  bcftools view -v snps -Oz -o "$out_chr_dir/snps_chr${chr}.vcf.gz" "$in_vcf"
  bcftools index -t "$out_chr_dir/snps_chr${chr}.vcf.gz"

  echo "Done chr${chr}: $out_chr_dir/snps_chr${chr}.vcf.gz"
done

echo "✓ VCFs with only SNP in pangenome_dir/SNP_dir/chr*/"
