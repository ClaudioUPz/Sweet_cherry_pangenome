#!/usr/bin/env bash
#SBATCH --job-name=snps_gwas
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16
#SBATCH --output=snps_gwas.%N.%j.out
#SBATCH --error=snps_gwas.%N.%j.err

set -euo pipefail

########################################################################
# Load modules
########################################################################
module load bcftools
module load samtools

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate rGWAS

########################################################################
# Parameters
########################################################################
DST="/path/to/project/8.GWAS/ref_and_bam_bai"
REF="${DST}/combined_reference.fa"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
KEEP_VCF="${KEEP_VCF:-0}"
LIST="${1:-}"

########################################################################
# Input checks
########################################################################
[[ -s "$REF" ]] || { echo "ERROR: missing reference FASTA: $REF"; exit 2; }
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

cd "$DST"

########################################################################
# Build BAM list
########################################################################
if [[ -n "$LIST" ]]; then
  cp -f "$LIST" bamlist.txt
else
  ls -1 *.merged.bam > bamlist.txt
fi

[[ -s bamlist.txt ]] || { echo "ERROR: no merged BAM files found in $DST"; exit 3; }

while read -r B; do
  [[ -s "${B}.bai" ]] || samtools index -@ "$THREADS" "$B"
done < bamlist.txt

########################################################################
# Output directories
########################################################################
JOBTAG="${SLURM_JOB_ID:-$$}"
OUTDIR="${DST}/mvp_snps"
TMPDIR="${DST}/snp_tmp_${JOBTAG}"
mkdir -p "$OUTDIR" "$TMPDIR"

FIL_VCF="${TMPDIR}/snps.cohort.filtered.vcf.gz"
MAT_TXT="${OUTDIR}/snp_matrix_num.txt"
MAP_TXT="${OUTDIR}/snp_map.txt"
MAP_NUM="${OUTDIR}/snp_map.numeric.txt"
IDMAP="${OUTDIR}/id_mapping_snps.csv"

########################################################################
# Cohort SNP calling and filtering
########################################################################
echo "[INFO] Running cohort SNP calling..."

bcftools mpileup -f "$REF" -b bamlist.txt \
    -q 20 -Q 20 -I -a DP,AD,SP \
    -Ou --threads "$THREADS" \
| bcftools call -m -v -Ou --threads "$THREADS" \
| bcftools norm -f "$REF" -m -both -Ou --threads "$THREADS" \
| bcftools view -v snps -Ou --threads "$THREADS" \
| bcftools sort -T "$TMPDIR" -m 4G -Ou \
| bcftools +fill-tags -Ou -- -t AN,AC,AF,F_MISSING \
| bcftools view -i 'QUAL>=30 && F_MISSING<0.2 && AF>=0.01 && AF<=0.99' \
    -Oz --threads "$THREADS" > "$FIL_VCF"

bcftools index -t -f "$FIL_VCF"

########################################################################
# Export rMVP matrix
########################################################################
echo "[INFO] Exporting rMVP files..."

{
  printf "ID"
  bcftools query -l "$FIL_VCF" | awk '{printf "\t"$0}'
  printf "\n"
} > "$MAT_TXT"

bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' "$FIL_VCF" \
| awk 'BEGIN{OFS="\t"}{
    chr=$1; pos=$2; id=$3; if(id=="."||id=="") id=chr":"pos;
    printf("%s", id);
    for(i=4;i<=NF;i++){
      g=$i; gsub(/\|/,"/",g);
      if(g=="0/0") v=0;
      else if(g=="0/1"||g=="1/0") v=1;
      else if(g=="1/1") v=2;
      else v="NA";
      printf("\t%s", v);
    }
    printf("\n");
  }' >> "$MAT_TXT"

echo -e "SNP\tChrom\tPosition" > "$MAP_TXT"
bcftools query -f '%CHROM\t%POS\t%ID\n' "$FIL_VCF" \
| awk 'BEGIN{OFS="\t"}{
    chr=$1; pos=$2; id=$3; if(id=="."||id=="") id=chr":"pos;
    print id, chr, pos
  }' >> "$MAP_TXT"

########################################################################
# Generate unique marker IDs and ID mapping table
########################################################################
echo -e "SNP_new,SNP_old,Chrom,Position" > "$IDMAP"

awk 'BEGIN{FS=OFS="\t"}
     NR==1{next}
     {
       old=$1; chr=$2; pos=$3;
       count[old]++;
       new=(count[old]==1 ? old : old"__dup"count[old]);
       print new","old","chr","pos >> "'"$IDMAP"'";
       rows[++n]=new OFS chr OFS pos;
     }
     END{
       print "SNP\tChrom\tPosition" > "'"$MAP_TXT"'.tmp";
       for(i=1;i<=n;i++) print rows[i] >> "'"$MAP_TXT"'.tmp";
     }' "$MAP_TXT"

mv -f "${MAP_TXT}.tmp" "$MAP_TXT"

########################################################################
# Generate numeric map
########################################################################
echo -e "SNP\tChrom\tPosition" > "$MAP_NUM"

awk 'BEGIN{FS=OFS="\t"}
     NR==1{next}
     {
       snp=$1; chr=$2; pos=$3;
       if (chr ~ /^chr[0-9]+__/) {
         sub(/^chr/,"",chr)
         sub(/__.*/,"",chr)
         c=chr+0
       } else {
         c=0
       }
       print snp, c, pos
     }' "$MAP_TXT" >> "$MAP_NUM"

########################################################################
# Optional cleanup
########################################################################
if [[ "${KEEP_VCF}" -eq 0 ]]; then
  echo "[INFO] Removing filtered VCF to save disk space"
  rm -f "$FIL_VCF" "${FIL_VCF}.tbi"
else
  echo "[INFO] Keeping filtered VCF: $FIL_VCF"
fi

find "$TMPDIR" -type f -name 'bcftools.*' -delete || true

########################################################################
# Summary
########################################################################
if [[ -s "$FIL_VCF" ]]; then
  nvars=$(bcftools view -H "$FIL_VCF" | wc -l | awk '{print $1}')
else
  nvars="NA"
fi

echo "[OK] SNP variants     : $nvars"
echo "[OK] rMVP matrix      : $MAT_TXT"
echo "[OK] rMVP map         : $MAP_TXT"
echo "[OK] rMVP numeric map : $MAP_NUM"
echo "[OK] ID mapping       : $IDMAP"
echo "[OK] Temporary files  : $TMPDIR"
