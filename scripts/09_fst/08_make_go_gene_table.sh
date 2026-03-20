#!/usr/bin/env bash
set -euo pipefail

GO_TSV="${1:-}"
GO2GENE_COORDS="${2:-}"
OUTFILE="${3:-}"

if [[ -z "$GO_TSV" || -z "$GO2GENE_COORDS" || -z "$OUTFILE" ]]; then
  echo "Usage:"
  echo "  bash $0 <go_enrichment.tsv> <go_to_gene.coords.tsv> <output_table.tsv>"
  exit 1
fi

[[ -f "$GO_TSV" ]] || { echo "ERROR: file not found: $GO_TSV"; exit 1; }
[[ -f "$GO2GENE_COORDS" ]] || { echo "ERROR: file not found: $GO2GENE_COORDS"; exit 1; }

for tool in awk mktemp rm; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing required tool '$tool'"; exit 1; }
done

tmp_go_desc="$(mktemp)"
trap 'rm -f "$tmp_go_desc"' EXIT

echo "Generating GO-gene summary table..."

awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="ID") id_col=i
      if($i=="Description") desc_col=i
    }
    if(!id_col || !desc_col){
      print "ERROR: required columns ID and/or Description not found in GO enrichment table header" > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    print $id_col, $desc_col
  }' "$GO_TSV" > "$tmp_go_desc"

awk 'BEGIN{FS=OFS="\t"}
  NR==FNR{
    desc[$1]=$2
    next
  }
  {
    go=$1
    gene=$2
    chr=$3
    start=$4
    end=$5
    strand=$6
    go_desc=(go in desc ? desc[go] : "NA")
    print go, go_desc, gene, chr, start, end, strand
  }' "$tmp_go_desc" "$GO2GENE_COORDS" > "$OUTFILE"

echo "Done"
echo "Output:"
echo "$OUTFILE"
