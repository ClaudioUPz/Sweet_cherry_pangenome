#!/usr/bin/env bash
set -euo pipefail

FST_OUTDIR="${1:-}"
GFF3="${2:-}"
OUTDIR="${3:-}"
MERGE_GAP="${4:-}"
GLM_SIGNALS_CSV="${5:-}"
FLANK_BP="${6:-}"
INTERPRO_TSV="${7:-}"
GWAS_SV_BED="${8:-}"

if [[ -z "$FST_OUTDIR" || -z "$GFF3" || -z "$OUTDIR" || -z "$MERGE_GAP" || -z "$GLM_SIGNALS_CSV" || -z "$FLANK_BP" || -z "$INTERPRO_TSV" ]]; then
  echo "Usage:"
  echo "  bash $0 <fst_outliers_dir> <annotation.gff3> <output_dir> <merge_gap_bp> <gwas_signals.csv> <flank_bp> <interpro.tsv> [gwas_sv_hits.bed]"
  exit 1
fi

mkdir -p "$OUTDIR"

for tool in bedtools awk sort cut tail grep comm wc cat; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing required tool '$tool'"; exit 1; }
done

[[ -d "$FST_OUTDIR" ]] || { echo "ERROR: directory not found: $FST_OUTDIR"; exit 1; }
[[ -f "$GFF3" ]] || { echo "ERROR: file not found: $GFF3"; exit 1; }
[[ -f "$GLM_SIGNALS_CSV" ]] || { echo "ERROR: file not found: $GLM_SIGNALS_CSV"; exit 1; }
[[ -f "$INTERPRO_TSV" ]] || { echo "ERROR: file not found: $INTERPRO_TSV"; exit 1; }

norm_chr_awk='function norm(c){ sub(/__.*/,"",c); return c }'

echo "[INFO] Output directory: $OUTDIR"
echo

FST_ALL_BED="$OUTDIR/fst_outliers_all.bed"
FST_MERGED_BED="$OUTDIR/fst_outliers_merged_d${MERGE_GAP}.bed"
FST_MERGED_CHROMFIX="$OUTDIR/fst_outliers_merged_d${MERGE_GAP}.chromfix.bed"

tmp_cat="$OUTDIR/.tmp_fst_cat.tsv"
> "$tmp_cat"

shopt -s nullglob
files=( "$FST_OUTDIR"/Fst_outliers_chr*.tsv )
shopt -u nullglob

if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: no files matching $FST_OUTDIR/Fst_outliers_chr*.tsv were found"
  exit 1
fi

for f in "${files[@]}"; do
  tail -n +2 "$f" >> "$tmp_cat"
done

awk 'BEGIN{OFS="\t"} NF>=3 {print $1, $2-1, $3}' "$tmp_cat" \
  | sort -k1,1 -k2,2n > "$FST_ALL_BED"

bedtools merge -i "$FST_ALL_BED" -d "$MERGE_GAP" > "$FST_MERGED_BED"

awk 'BEGIN{OFS="\t"} '"$norm_chr_awk"' {print norm($1), $2, $3}' "$FST_MERGED_BED" \
  | sort -k1,1 -k2,2n > "$FST_MERGED_CHROMFIX"

echo "[OK] FST merged BED: $FST_MERGED_CHROMFIX"
echo

GENE_BED="$OUTDIR/genes.bed"
ALL_GENES_LIST="$OUTDIR/all_genes.list.txt"

awk 'BEGIN{FS="\t"; OFS="\t"}
  $0 !~ /^#/ && $3=="gene" {
    chr=$1; start=$4-1; end=$5; strand=$7;
    id="NA";
    n=split($9,a,";");
    for(i=1;i<=n;i++){
      if(a[i] ~ /^ID=/){sub(/^ID=/,"",a[i]); id=a[i]}
      else if(id=="NA" && a[i] ~ /^Name=/){sub(/^Name=/,"",a[i]); id=a[i]}
    }
    print chr, start, end, id, strand
  }' "$GFF3" \
  | sort -k1,1 -k2,2n > "$GENE_BED"

cut -f4 "$GENE_BED" | sort -u > "$ALL_GENES_LIST"

echo "[OK] Total genes in annotation: $(wc -l < "$ALL_GENES_LIST")"
echo

GENES_GLOBAL_LIST="$OUTDIR/genes_in_fst_outliers.list.txt"

bedtools intersect -a "$GENE_BED" -b "$FST_MERGED_CHROMFIX" -wa -u \
  | cut -f4 | sort -u > "$GENES_GLOBAL_LIST"

echo "[OK] Genes in FST outlier regions: $(wc -l < "$GENES_GLOBAL_LIST")"
echo

GWAS_SNP_HITS="$OUTDIR/gwas_snp_hits.tsv"

awk -F',' 'BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      x=$i; gsub(/"/,"",x);
      if(x=="Chrom") c=i;
      if(x=="Position") p=i;
      if(x ~ /Trait\.GLM/ || x=="P" || x=="pvalue" || x=="p.value") pv=i;
    }
    if(!c || !p){
      print "ERROR: required columns Chrom and/or Position were not found in GWAS header" > "/dev/stderr";
      exit 1
    }
    print "CHROM","POS","P";
    next
  }
  {
    chrom=$c; pos=$p;
    gsub(/"/,"",chrom);
    gsub(/"/,"",pos);
    pval="NA";
    if(pv){
      pval=$(pv);
      gsub(/"/,"",pval)
    }
    if(chrom != "" && pos != "") print chrom, pos, pval;
  }' "$GLM_SIGNALS_CSV" > "$GWAS_SNP_HITS"

echo "[OK] GWAS SNP hits: $(($(wc -l < "$GWAS_SNP_HITS") - 1))"
echo

GWAS_SNP_BED="$OUTDIR/gwas_snp_pm${FLANK_BP}.chromfix.merged.bed"
GWAS_UNION_BED="$OUTDIR/gwas_union_pm${FLANK_BP}.chromfix.merged.bed"
FST_X_GWAS_BED="$OUTDIR/fst_outliers_x_gwas.chromfix.bed"
GENES_FOCAL_LIST="$OUTDIR/genes_in_fst_x_gwas.list.txt"

awk 'BEGIN{OFS="\t"} '"$norm_chr_awk"' NR==1{next}
  {
    c=norm($1);
    pos=$2;
    start=pos-FL;
    if(start<1) start=1;
    end=pos+FL;
    print c, start-1, end
  }' FL="$FLANK_BP" "$GWAS_SNP_HITS" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$GWAS_SNP_BED"

tmp_union="$OUTDIR/.tmp_gwas_union.bed"
cat "$GWAS_SNP_BED" > "$tmp_union"

if [[ -n "${GWAS_SV_BED:-}" ]]; then
  if [[ -f "$GWAS_SV_BED" ]]; then
    GWAS_SV_BED_OUT="$OUTDIR/gwas_sv_pm${FLANK_BP}.chromfix.merged.bed"
    awk 'BEGIN{OFS="\t"} '"$norm_chr_awk"' {
      c=norm($1);
      s=$2-FL; if(s<0) s=0;
      e=$3+FL;
      print c, s, e
    }' FL="$FLANK_BP" "$GWAS_SV_BED" \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - > "$GWAS_SV_BED_OUT"
    cat "$GWAS_SV_BED_OUT" >> "$tmp_union"
  else
    echo "[WARN] Optional GWAS SV BED file not found: $GWAS_SV_BED"
  fi
fi

sort -k1,1 -k2,2n "$tmp_union" | bedtools merge -i - > "$GWAS_UNION_BED"

bedtools intersect -a "$FST_MERGED_CHROMFIX" -b "$GWAS_UNION_BED" > "$FST_X_GWAS_BED"

bedtools intersect -a "$GENE_BED" -b "$FST_X_GWAS_BED" -wa -u \
  | cut -f4 | sort -u > "$GENES_FOCAL_LIST"

echo "[OK] Genes in FST-GWAS overlapping regions: $(wc -l < "$GENES_FOCAL_LIST")"
echo

G2GO_TX="$OUTDIR/gene2go.transcript.tsv"
G2GO_GENE="$OUTDIR/gene2go.gene.tsv"

awk 'BEGIN{FS="\t"; OFS="\t"}
{
  id=$1;
  line=$0;
  gos="";
  while(match(line, /GO:[0-9]{7}/)){
    term=substr(line, RSTART, RLENGTH);
    if(gos=="") gos=term; else gos=gos";"term;
    line=substr(line, RSTART+RLENGTH);
  }
  if(gos!="") print id, gos;
}' "$INTERPRO_TSV" \
| awk 'BEGIN{FS="\t"; OFS="\t"}
  {
    id=$1; gos=$2;
    n=split(gos,a,";");
    for(i=1;i<=n;i++){
      if(a[i]!="") seen[id SUBSEP a[i]]=1;
    }
  }
  END{
    for(k in seen){
      split(k,parts,SUBSEP);
      id=parts[1];
      go=parts[2];
      if(id in out) out[id]=out[id]";"go; else out[id]=go;
    }
    for(id in out) print id, out[id];
  }' \
| sort -u > "$G2GO_TX"

awk 'BEGIN{FS="\t"; OFS="\t"}
  {
    tx=$1; gos=$2;
    gene=tx;
    sub(/\.[0-9]+$/,"",gene);
    n=split(gos,a,";");
    for(i=1;i<=n;i++){
      if(a[i]!="") seen[gene SUBSEP a[i]]=1;
    }
  }
  END{
    for(k in seen){
      split(k,parts,SUBSEP);
      g=parts[1];
      go=parts[2];
      if(g in out) out[g]=out[g]";"go; else out[g]=go;
    }
    for(g in out) print g, out[g];
  }' "$G2GO_TX" \
| sort -u > "$G2GO_GENE"

echo "[OK] Transcript-level gene2go entries: $(wc -l < "$G2GO_TX")"
echo "[OK] Gene-level gene2go entries: $(wc -l < "$G2GO_GENE")"
echo

GLOBAL_WITH_GO=$(comm -12 <(sort -u "$GENES_GLOBAL_LIST") <(cut -f1 "$G2GO_GENE" | sort -u) | wc -l)
FOCAL_WITH_GO=$(comm -12 <(sort -u "$GENES_FOCAL_LIST") <(cut -f1 "$G2GO_GENE" | sort -u) | wc -l)

echo "[INFO] Genes in FST outlier regions with GO annotations: $GLOBAL_WITH_GO / $(wc -l < "$GENES_GLOBAL_LIST")"
echo "[INFO] Genes in FST-GWAS overlapping regions with GO annotations: $FOCAL_WITH_GO / $(wc -l < "$GENES_FOCAL_LIST")"
echo
echo "Done."
echo "gene2go file: $G2GO_GENE"
echo "FST outlier genes: $GENES_GLOBAL_LIST"
echo "FST-GWAS genes: $GENES_FOCAL_LIST"
echo "Background gene list: $ALL_GENES_LIST"
