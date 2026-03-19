#!/bin/bash
#SBATCH --job-name=validate_d1_vg
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/validate_d1_vg.%j.out
#SBATCH --error=logs/validate_d1_vg.%j.err

module load vg

echo "Starting validation of .d1.vg files"
date

BASE_DIR="/path/to/workdir"

mkdir -p logs

for i in {1..8}; do
  CHR="chr${i}"
  VG_DIR="${BASE_DIR}/${CHR}/output/${CHR}_run.chroms"

  echo "Processing: $VG_DIR"
  if [ -d "$VG_DIR" ]; then
    for vg_file in "$VG_DIR"/*.d1.vg; do
      if [ -f "$vg_file" ]; then
        echo "Validating $(basename "$vg_file")"
        vg validate "$vg_file"
      else
        echo "No .d1.vg files found in $VG_DIR"
      fi
    done
  else
    echo "Directory not found: $VG_DIR"
  fi
  echo "--------------------------------------"
done

echo "Validation completed"
date
