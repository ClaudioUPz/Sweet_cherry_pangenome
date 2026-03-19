#!/bin/bash
#SBATCH --job-name=snps_vcf
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -o logs/snp_vcf.%j.out
#SBATCH -e logs/snp_vcf.%j.err

set -euo pipefail

# Cargar módulo
module load bcftools

BASE="/home/curra/pangenome"
IN_DIR="$BASE/5.pangenome"
OUT_DIR="$BASE/6.SNP"

mkdir -p "$OUT_DIR" logs

for chr in {1..8}; do
  in_vcf="$IN_DIR/chr${chr}/output/chr${chr}_run.vcf.gz"
  if [[ ! -f "$in_vcf" ]]; then
    echo "ERROR: No existe $in_vcf"; exit 1
  fi

  out_chr_dir="$OUT_DIR/chr${chr}"
  mkdir -p "$out_chr_dir"

  # SNPs (mantiene bi y multi-alélicos). Indexa con tabix.
  bcftools view -v snps -Oz -o "$out_chr_dir/snps_chr${chr}.vcf.gz" "$in_vcf"
  bcftools index -t "$out_chr_dir/snps_chr${chr}.vcf.gz"

  # Si SOLO dejamos bialélicos, descomentar la línea siguiente en lugar de la anterior:
  # bcftools view -v snps -m2 -M2 -Oz -o "$out_chr_dir/snps_chr${chr}.vcf.gz" "$in_vcf" && bcftools index -t "$out_chr_dir/snps_chr${chr}.vcf.gz"

  echo "Listo chr${chr}: $out_chr_dir/snps_chr${chr}.vcf.gz"
done

echo "✓ VCFs con solo SNPs en pangenome/6.SNP/chr*/"
