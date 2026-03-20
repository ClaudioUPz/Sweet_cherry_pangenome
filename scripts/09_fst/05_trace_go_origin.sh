#!/usr/bin/env bash
set -euo pipefail

GO_TSV="${1:-}"
GENES_BED="${2:-}"
FST_OUTLIERS_BED="${3:-}"
OUTDIR="${4:-}"

if [[ -z "$GO_TSV" || -z "$GENES_BED" || -z "$FST_OUTLIERS_BED" || -z "$OUTDIR" ]]; then
  echo "Usage:"
  echo "  bash $0 <go_enrichment.tsv> <genes.bed> <fst_outliers_merged.bed> <output_dir>"
  exit 1
fi

[[ -f "$GO_TSV" ]] || { echo "ERROR: file not found: $GO_TSV"; exit 1; }
[[ -f "$GENES_BED" ]] || { echo "ERROR: file not found: $GENES_BED"; exit 1; }
[[ -f "$FST_OUTLIERS_BED" ]] || { echo "ERROR: file not found: $FST_OUTLIERS_BED"; exit 1; }

for tool in bedtools awk sort uniq cut head cat wc; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing required tool '$tool'"; exit 1; }
done

mkdir -p "$OUTDIR"

GO2GENE="$OUTDIR/go_to_gene.tsv"
GO2GENE_COORDS="$OUTDIR/go_to_gene.coords.tsv"
BY_CHR="$OUTDIR/go_by_chr.tsv"
SPAN_BY_CHR="$OUTDIR/go_span_by_chr.tsv"
GENES_IN_OUTLIERS="$OUTDIR/go_genes_in_outlier_windows.tsv"
GO_WINDOWS="$OUTDIR/go_windows.tsv"
SUMMARY="$OUTDIR/go_summary.tsv"

echo "[INFO] GO enrichment table : $GO_TSV"
echo "[INFO] Gene coordinates    : $GENES_BED"
echo "[INFO] FST outlier BED     : $FST_OUTLIERS_BED"
echo "[INFO] Output directory    : $OUTDIR"
echo

echo "[1/7] Expanding GO enrichment table to GO-to-gene mapping..."

awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="ID") go_col=i
      if($i=="geneID") gene_col=i
    }
    if(!go_col || !gene_col){
      print "ERROR: required columns ID and/or geneID not found in GO table header" > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    go=$go_col
    genes=$gene_col
    if(genes=="" || genes=="NA") next
    n=split(genes,a,"/")
    for(i=1;i<=n;i++){
      if(a[i]!="") print go,a[i]
    }
  }' "$GO_TSV" \
| sort -u > "$GO2GENE"

echo "[OK] GO-to-gene table: $GO2GENE (n=$(wc -l < "$GO2GENE"))"
echo

echo "[2/7] Mapping GO-associated genes to genomic coordinates..."

awk 'BEGIN{FS="\t"; OFS="\t"}
  NR==FNR {coord[$4]=$1"\t"$2"\t"$3"\t"$5; next}
  {
    go=$1
    gene=$2
    if(gene in coord) print go,gene,coord[gene]
    else print go,gene,"NA","NA","NA","NA"
  }' "$GENES_BED" "$GO2GENE" \
> "$GO2GENE_COORDS"

echo "[OK] GO gene coordinates: $GO2GENE_COORDS"
echo

echo "[3/7] Summarizing GO-associated genes by chromosome..."

awk 'BEGIN{FS="\t"; OFS="\t"}
  $3!="NA" {count[$1 SUBSEP $3]++}
  END{
    print "GO","CHR","N_genes"
    for(k in count){
      split(k,a,SUBSEP)
      print a[1],a[2],count[k]
    }
  }' "$GO2GENE_COORDS" \
| sort -k1,1 -k3,3nr > "$BY_CHR"

echo "[OK] GO by chromosome: $BY_CHR"
echo

echo "[4/7] Calculating genomic span of GO-associated genes by chromosome..."

awk 'BEGIN{FS="\t"; OFS="\t"}
  $3!="NA"{
    go=$1
    chr=$3
    s=$4
    e=$5
    key=go SUBSEP chr
    if(!(key in minS) || s<minS[key]) minS[key]=s
    if(!(key in maxE) || e>maxE[key]) maxE[key]=e
    n[key]++
  }
  END{
    print "GO","CHR","N_genes","min_start0","max_end","span_bp"
    for(k in n){
      split(k,a,SUBSEP)
      span=maxE[k]-minS[k]
      print a[1],a[2],n[k],minS[k],maxE[k],span
    }
  }' "$GO2GENE_COORDS" \
| sort -k1,1 -k3,3nr > "$SPAN_BY_CHR"

echo "[OK] GO genomic span summary: $SPAN_BY_CHR"
echo

echo "[5/7] Intersecting GO-associated genes with FST outlier windows..."

TMP_GO_GENES_BED="$OUTDIR/.tmp_go_genes_coords.bed"

awk 'BEGIN{FS="\t"; OFS="\t"}
  $3!="NA" {print $3,$4,$5,$2,$1}' "$GO2GENE_COORDS" \
| sort -k1,1 -k2,2n > "$TMP_GO_GENES_BED"

bedtools intersect -wa -wb \
  -a "$TMP_GO_GENES_BED" \
  -b "$FST_OUTLIERS_BED" \
> "$GENES_IN_OUTLIERS"

echo "[OK] GO genes intersecting FST outlier windows: $GENES_IN_OUTLIERS (n=$(wc -l < "$GENES_IN_OUTLIERS"))"
echo

echo "[6/7] Summarizing unique FST outlier windows contributing to each GO term..."

awk 'BEGIN{FS="\t"; OFS="\t"}
  {
    chr=$1
    wstart=$6
    wend=$7
    go=$5
    win=chr":"wstart"-"wend
    seen[go SUBSEP win]=1
    chrseen[go SUBSEP chr]=1
    geneset[go SUBSEP $4]=1
  }
  END{
    print "GO","N_unique_windows","CHR_list","N_unique_genes_in_outliers"
    for(k in seen){
      split(k,a,SUBSEP)
      go=a[1]
      nwin[go]++
    }
    for(k in geneset){
      split(k,a,SUBSEP)
      go=a[1]
      ngenes[go]++
    }
    for(k in chrseen){
      split(k,a,SUBSEP)
      go=a[1]
      chr=a[2]
      if(go in chrs) chrs[go]=chrs[go]","chr
      else chrs[go]=chr
    }
    for(go in nwin){
      print go,nwin[go],chrs[go],ngenes[go]
    }
  }' "$GENES_IN_OUTLIERS" \
| sort -k2,2nr > "$GO_WINDOWS"

echo "[OK] GO window summary: $GO_WINDOWS"
echo

echo "[7/7] Building final GO summary table..."

TMP_GO_STATS="$OUTDIR/.tmp_go_stats.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="ID") go_col=i
      if($i=="p.adjust") padj_col=i
      if($i=="Count") count_col=i
    }
    if(!go_col || !padj_col || !count_col){
      print "ERROR: required columns ID, p.adjust, and/or Count not found in GO table header" > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    print $go_col,$padj_col,$count_col
  }' "$GO_TSV" > "$TMP_GO_STATS"

TMP_GO_NGENES="$OUTDIR/.tmp_go_ngenes.tsv"

awk 'BEGIN{FS="\t"; OFS="\t"}
  $3!="NA" {seen[$1 SUBSEP $2]=1}
  END{
    for(k in seen){
      split(k,a,SUBSEP)
      go=a[1]
      n[go]++
    }
    for(go in n) print go,n[go]
  }' "$GO2GENE_COORDS" > "$TMP_GO_NGENES"

awk 'BEGIN{FS="\t"; OFS="\t"}
  NR==FNR {ng[$1]=$2; next}
  {print $1,$2,$3,(($1 in ng)?ng[$1]:"NA")}' "$TMP_GO_NGENES" "$TMP_GO_STATS" \
> "$OUTDIR/.tmp_go_stats_ng.tsv"

awk 'BEGIN{FS="\t"; OFS="\t"}
  NR==FNR {win[$1]=$2; chrs[$1]=$3; ngenout[$1]=$4; next}
  {
    go=$1
    print $0,((go in win)?win[go]:"0"),((go in chrs)?chrs[go]:"NA"),((go in ngenout)?ngenout[go]:"0")
  }' "$GO_WINDOWS" "$OUTDIR/.tmp_go_stats_ng.tsv" \
| awk 'BEGIN{FS="\t"; OFS="\t"}
  NR==1{
    print "GO","p_adjust","Count","N_genes_total","N_unique_windows_in_FST_outliers","CHR_list_in_FST_outliers","N_unique_genes_in_FST_outliers"
  }
  {print $1,$2,$3,$4,$5,$6,$7}' \
> "$SUMMARY"

echo "[OK] Final GO summary: $SUMMARY"
echo
echo "Done."
echo
echo "Preview of final GO summary:"
head -n 20 "$SUMMARY"
echo
echo "Preview of GO contribution by chromosome:"
head -n 40 "$BY_CHR"
echo
echo "Preview of GO genomic span summary:"
head -n 40 "$SPAN_BY_CHR"
