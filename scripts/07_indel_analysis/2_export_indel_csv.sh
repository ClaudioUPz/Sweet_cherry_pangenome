#!/bin/bash
#SBATCH --job-name=indels_csv
#SBATCH --time=02:00:00
#SBATCH -o indels_csv.%N.%j.out
#SBATCH -e indels_csv.%N.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

# --- Config ---
BASE="/path/to/project"
CHRS="${CHRS:-1 2 3 4 5 6 7 8}"

# --- TOOLS ---
module load bcftools
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found"; exit 1; }

for chr in ${CHRS}; do
  vcf="${BASE}/7.INDEL/chr${chr}/indels_chr${chr}.vcf.gz"
  out="${BASE}/7.INDEL/chr${chr}/indels_details_chr${chr}.csv"

  [[ -s "${vcf}" ]] || { echo "ERROR: it is absent ${vcf} (run Step 1)"; exit 1; }

  echo ">> chr${chr} -> $(basename "${out}")"
  echo "CHROM,POS,REF,ALT,SAMPLE,GT" > "${out}"

  # Line by (variant, sample) with ALT (GT contains 1..9)
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%GT]\n' "${vcf}" \
  | awk -F'\t' '{
      CH=$1; POS=$2; REF=$3; ALT=$4;
      for(i=5;i<=NF;i+=2){
        S=$(i); GT=$(i+1);
        if(GT ~ /[1-9]/) printf "%s,%s,%s,%s,%s,%s\n", CH,POS,REF,ALT,S,GT;
      }
    }' >> "${out}"

  nl=$(wc -l < "${out}" || echo 0)
  echo "OK chr${chr}: ${nl} lines (including header)"
done

echo "Step 2 done: 7.INDEL/chrN/indels_details_chrN.csv"
