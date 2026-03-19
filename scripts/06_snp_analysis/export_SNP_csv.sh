#!/bin/bash
#SBATCH --job-name=snps_csv
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -o logs/snp_csv.%j.out
#SBATCH -e logs/snp_csv.%j.err

set -euo pipefail
module load bcftools

BASE="/home/curra/pangenome"
SNP_DIR="$BASE/6.SNP"

mkdir -p logs

for chr in {1..8}; do
  vcf="$SNP_DIR/chr${chr}/snps_chr${chr}.vcf.gz"
  [[ -f "$vcf" ]] || { echo "Falta $vcf (corre 01 primero)"; exit 1; }

  outdir="$SNP_DIR/chr${chr}"
  names_file="$outdir/samples.txt"
  bcftools query -l "$vcf" > "$names_file"

  tmp="$outdir/_gt.tsv"
  # Matriz: CHROM POS REF ALT [GT por muestra...]
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$vcf" > "$tmp"

  # CSV largo: SOLO filas donde el GT de la muestra contiene alelo ALT (1–9)
  awk -v OFS=',' -v names="$names_file" '
    BEGIN{
      i=1; while((getline n < names)>0){ sample[i]=n; i++ } close(names);
      print "CHROM","POS","REF","ALT","SAMPLE","GT";
    }
    {
      for(j=5;j<=NF;j++){
        gt=$j
        if (gt ~ /[1-9]/) {           # incluye 0/1, 1/1, 1|2, 1, 2/2, etc.; excluye 0, 0/0, ., ./., .|.
          idx=j-4; s=sample[idx];
          print $1,$2,$3,$4,s,gt;     # CHROM queda tal cual (p.ej. chr1__ReginaC)
        }
      }
    }' "$tmp" > "$outdir/snps_details_chr${chr}.csv"

  rm -f "$tmp"
  echo "✓ $outdir/snps_details_chr${chr}.csv"
done

echo "Listo."
