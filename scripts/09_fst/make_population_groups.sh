#!/usr/bin/env bash

set -euo pipefail

INPUT_TSV="${1:-sample_groups_checked.tsv}"
OUTPUT_TSV="${2:-population_groups.tsv}"

echo -e "#sample\tpopulation" > "$OUTPUT_TSV"

tail -n +2 "$INPUT_TSV" | while read -r SAMPLE STATUS; do
  case "$STATUS" in
    landrace)
      echo -e "${SAMPLE}\tlandrace" >> "$OUTPUT_TSV"
      ;;
    modern_breeding)
      echo -e "${SAMPLE}\tmodern" >> "$OUTPUT_TSV"
      ;;
    early_selection)
      echo -e "${SAMPLE}\tearly" >> "$OUTPUT_TSV"
      ;;
    *)
      echo "WARNING: unknown group for sample $SAMPLE: $STATUS" >&2
      ;;
  esac
done

echo "Population group table written to: $OUTPUT_TSV"
