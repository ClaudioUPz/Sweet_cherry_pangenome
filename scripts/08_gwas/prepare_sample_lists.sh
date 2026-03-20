#!/usr/bin/env bash

# ============================================================
# Script: prepare_sample_lists.sh
# Description:
# Generate ordered lists of sample IDs and paired-end FASTQ files
# for downstream read mapping (VG / Minimap2).
#
# Outputs:
# - ID_list.txt
# - R1_list.txt
# - R2_list.txt
#
# Requirements:
# - Directory structure: one folder per sample (e.g., VXXXX)
# - FASTQ files containing "_1" (R1) and "_2" (R2)
# ============================================================

set -euo pipefail

BASE_DIR="/path/to/project/2.short_reads"
cd "$BASE_DIR" || { echo "ERROR: directory not found: $BASE_DIR"; exit 1; }

# Output files
ID_LIST="ID_list.txt"
R1_LIST="R1_list.txt"
R2_LIST="R2_list.txt"

# Initialize output files
> "$ID_LIST"
> "$R1_LIST"
> "$R2_LIST"

echo "[INFO] Scanning sample directories..."

# Get sorted list of sample directories (e.g., VXXXX)
find . -maxdepth 1 -type d -regex './V[0-9]+' -printf '%f\n' | sort > "$ID_LIST"

echo "[INFO] Generating FASTQ lists..."

# Iterate through samples
while read -r sample; do

    # R1 file
    R1=$(ls "$sample"/*_1*.fq.gz 2>/dev/null | head -n 1 || true)
    if [[ -n "$R1" ]]; then
        realpath "$R1" >> "$R1_LIST"
    else
        echo "[WARNING] Missing R1 file for $sample" >&2
        echo "" >> "$R1_LIST"
    fi

    # R2 file
    R2=$(ls "$sample"/*_2*.fq.gz 2>/dev/null | head -n 1 || true)
    if [[ -n "$R2" ]]; then
        realpath "$R2" >> "$R2_LIST"
    else
        echo "[WARNING] Missing R2 file for $sample" >&2
        echo "" >> "$R2_LIST"
    fi

done < "$ID_LIST"

echo "[INFO] Sample lists generated:"
echo " - $ID_LIST"
echo " - $R1_LIST"
echo " - $R2_LIST"
