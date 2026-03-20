#!/usr/bin/env bash
#SBATCH --job-name=sv_calling_delly
#SBATCH --mem=264G
#SBATCH --cpus-per-task=16
#SBATCH --output=sv_calling.%N.%j.out
#SBATCH --error=sv_calling.%N.%j.err

set -euo pipefail
trap 'echo "[FATAL] line $LINENO: $BASH_COMMAND (code=$?)" >&2' ERR

########################################################################
# Load modules
########################################################################
module load vg
module load delly
module load bcftools
module load samtools

CLEAN_GAM="${CLEAN_GAM:-1}"

########################################################################
# Input/output paths
########################################################################
SRC="/path/to/project/5.pangenome"
DST="/path/to/project/8.GWAS/ref_and_bam_bai"
LIST="${1:-bam_list.txt}"

mkdir -p "$DST"
mkdir -p "$DST/tmp"
export TMPDIR="$DST/tmp"

########################################################################
# Graph preparation
########################################################################
GBZ="${SRC}/merged_graph.gbz"
XG="${DST}/merged_graph.xg"

[[ -s "$GBZ" ]] || { echo "ERROR: missing GBZ graph"; exit 1; }

if [[ ! -s "$XG" ]]; then
  echo "Converting GBZ to XG..."
  vg convert -x "$GBZ" > "$XG"
fi

########################################################################
# Extract reference paths and build FASTA
########################################################################
mapfile -t REF_PATHS < <(
  vg paths -L -x "$XG" | grep -E '^REFERENCE#0#chr[1-8]__REFERENCE$' | sort
)

[[ ${#REF_PATHS[@]} -gt 0 ]] || { echo "ERROR: no reference paths found"; exit 1; }

SURJECT_P=()
for p in "${REF_PATHS[@]}"; do SURJECT_P+=(-p "$p"); done

PATHS_FILE="${DST}/reference_paths.txt"
printf "%s\n" "${REF_PATHS[@]}" > "$PATHS_FILE"

REF_FA="${DST}/reference.fa"
if [[ ! -s "$REF_FA" ]]; then
  vg paths -F -x "$XG" -p "$PATHS_FILE" | sed -E 's/^>REFERENCE#0#/>/' > "$REF_FA"
fi

[[ -s "${REF_FA}.fai" ]] || samtools faidx "$REF_FA"

########################################################################
# Functions
########################################################################
inject_if_needed() {
  local bam="$1"
  local name
  name=$(basename "$bam" .bam)

  local gam="${DST}/${name}.gam"

  if [[ ! -s "$gam" ]]; then
    vg inject -x "$XG" -t "$SLURM_CPUS_PER_TASK" "$bam" > "$gam"
  fi

  printf '%s\n' "$gam"
}

surject_one() {
  local bam="$1"
  local name
  name=$(basename "$bam" .bam)

  local sample="${name%__vggiraffe__combined_clean}"
  local root
  root="$(sed -E 's/-(1|2)$//' <<< "$sample")"

  local outbam="${DST}/${sample}.sorted.bam"

  if [[ -s "$outbam" && -s "${outbam}.bai" ]]; then
    printf '%s\t%s\n' "$root" "$outbam"
    return 0
  fi

  local gam
  gam="$(inject_if_needed "$bam")"

  vg surject -x "$XG" -b -t "$SLURM_CPUS_PER_TASK" "${SURJECT_P[@]}" "$gam" \
  | samtools view -h - \
  | awk 'BEGIN{FS=OFS="\t"}
         /^@SQ/ { sub(/\tSN:REFERENCE#0#/, "\tSN:"); print; next }
         /^@/   { print; next }
         { gsub(/^REFERENCE#0#/, "", $3); print }' \
  | samtools view -b - \
  | samtools addreplacerg -r "ID:${sample}" -r "SM:${root}" -r "PL:ILLUMINA" - \
  | samtools sort -T "$DST/tmp/${sample}" -o "$outbam" -

  samtools index "$outbam"

  if [[ "$CLEAN_GAM" -eq 1 ]]; then
    rm -f "$gam"
  fi

  printf '%s\t%s\n' "$root" "$outbam"
}

########################################################################
# Surjection and grouping
########################################################################
declare -A GROUP

while IFS= read -r bam; do
  [[ -s "$bam" ]] || continue

  if line="$(surject_one "$bam")"; then
    root="${line%%$'\t'*}"
    sbam="${line#*$'\t'}"
    GROUP["$root"]+="$sbam "
  fi
done < "$LIST"

########################################################################
# Merge replicates
########################################################################
MERGED_BAMS=()

for root in "${!GROUP[@]}"; do
  IFS=' ' read -r -a files <<< "${GROUP[$root]}"

  merged="${DST}/${root}.merged.bam"

  if (( ${#files[@]} == 1 )); then
    ln -sfn "${files[0]}" "$merged"
    ln -sfn "${files[0]}.bai" "$merged.bai"
  else
    samtools merge -@ "$SLURM_CPUS_PER_TASK" -o "$merged" "${files[@]}"
    samtools index "$merged"
  fi

  MERGED_BAMS+=("$merged")
done

########################################################################
# DELLY SV calling
########################################################################
BCFS=()

for mbam in "${MERGED_BAMS[@]}"; do
  root=$(basename "$mbam" .merged.bam)
  out="${DST}/${root}.delly.bcf"

  if [[ ! -s "$out" ]]; then
    delly call -g "$REF_FA" -o "$out" "$mbam"
  fi

  BCFS+=("$out")
done

########################################################################
# Merge SV sites
########################################################################
SITES="${DST}/delly.sites.bcf"
delly merge -o "$SITES" "${BCFS[@]}"

########################################################################
# Genotyping across samples
########################################################################
GENO_BCFS=()

for mbam in "${MERGED_BAMS[@]}"; do
  root=$(basename "$mbam" .merged.bam)
  out="${DST}/${root}.delly.geno.bcf"

  if [[ ! -s "$out" ]]; then
    delly call -g "$REF_FA" -v "$SITES" -o "$out" "$mbam"
  fi

  GENO_BCFS+=("$out")
done

########################################################################
# Final merged VCF
########################################################################
MERGED="${DST}/delly.merged.vcf.gz"

bcftools merge -m all -Oz -o "$MERGED" "${GENO_BCFS[@]}"
bcftools index -t "$MERGED"

echo "SV calling completed successfully"
