#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage:"
  echo "  bash make_group_lists.sh <sample_groups_checked.tsv> <output_dir>"
  exit 1
fi

SAMPLE_GROUPS=$1
OUTDIR=$2

echo "Sample group table: $SAMPLE_GROUPS"
echo "Output directory:   $OUTDIR"
echo

mkdir -p "$OUTDIR"

LANDRACE="$OUTDIR/landrace_samples.txt"
MODERN="$OUTDIR/modern_breeding_samples.txt"
EARLY="$OUTDIR/early_selection_samples.txt"

> "$LANDRACE"
> "$MODERN"
> "$EARLY"

echo "Generating sample lists by group..."

tail -n +2 "$SAMPLE_GROUPS" | while read -r SAMPLE STATUS; do
  case "$STATUS" in
    landrace)
      echo "$SAMPLE" >> "$LANDRACE"
      ;;
    modern_breeding)
      echo "$SAMPLE" >> "$MODERN"
      ;;
    early_selection)
      echo "$SAMPLE" >> "$EARLY"
      ;;
    *)
      echo "WARNING: unknown group for sample $SAMPLE: $STATUS" >&2
      ;;
  esac
done

echo
echo "Summary:"
echo "  Landrace:         $(wc -l < "$LANDRACE") samples"
echo "  Modern breeding:  $(wc -l < "$MODERN") samples"
echo "  Early selection:  $(wc -l < "$EARLY") samples"
echo
echo "Files written:"
echo "  $LANDRACE"
echo "  $MODERN"
echo "  $EARLY"
