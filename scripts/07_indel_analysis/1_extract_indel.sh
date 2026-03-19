#!/bin/bash
#SBATCH --job-name=indels_extract
#SBATCH --time=02:00:00
#SBATCH -o indels_extract.%N.%j.out
#SBATCH -e indels_extract.%N.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

# --- Config ---
BASE="/path/to/project"
# Choose chromosomes to process from environment variable CHRS (p.ej. CHRS="1 3 5")
CHRS="${CHRS:-1 2 3 4 5 6 7 8}"
# MODE it can be: auto (prefer final), final (obliges run.vcf.gz), raw (obliges run.raw.vcf.gz)
MODE="${MODE:-auto}"

# --- TOOLS ---
# Load bcftools
module load bcftools
command -v bcftools >/dev/null || { echo "ERROR: bcftools no encontrado en PATH ni módulos."; exit 1; }

echo "Usando MODE=${MODE}; cromosomas: ${CHRS}"
echo "BASE = ${BASE}"

for chr in ${CHRS}; do
  echo ">> chr${chr}"

  final_vcf="${BASE}/pangenome_dir/chr${chr}/output/chr${chr}_run.vcf.gz"
  raw_vcf="${BASE}/pangenome_dir/chr${chr}/output/chr${chr}_run.raw.vcf.gz"

  # Select the source VCF
  case "${MODE}" in
    final)
      SRC="${final_vcf}"
      [[ -s "${SRC}" ]] || { echo "ERROR: no existe ${SRC}"; exit 1; }
      ;;
    raw)
      SRC="${raw_vcf}"
      [[ -s "${SRC}" ]] || { echo "ERROR: no existe ${SRC}"; exit 1; }
      ;;
    auto|*)
      if [[ -s "${final_vcf}" ]]; then
        SRC="${final_vcf}"
      elif [[ -s "${raw_vcf}" ]]; then
        SRC="${raw_vcf}"
        echo "AVISO: using RAW for chr${chr} (didn't find final)."
      else
        echo "ERROR: Didn't find source VCF (final or raw) for chr${chr}."
        exit 1
      fi
      ;;
  esac

  outdir="${BASE}/7.INDEL/chr${chr}"
  mkdir -p "${outdir}"

  outvcf="${outdir}/indels_chr${chr}.vcf.gz"

  # Filter only INDELs and compress en BGZF

  # Index with bcftools (genera .tbi estilo tabix)
  bcftools index -f -t "${outvcf}"

  # Save the samples list
  bcftools query -l "${outvcf}" > "${outdir}/samples.txt"

  # Summary
  nvars=$(bcftools view -H "${outvcf}" | wc -l || echo 0)
  nsamp=$(wc -l < "${outdir}/samples.txt" || echo 0)
  echo "OK chr${chr}: ${nvars} INDELs, ${nsamp} muestras"
done

echo "Done: VCFs in ${BASE}/7.INDEL/chrN/indels_chrN.vcf.gz (+ .tbi) y samples.txt"
