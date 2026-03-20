#!/bin/bash
#SBATCH --job-name=vg_minimap2_compare
#SBATCH --time=240:00:00
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=20
#SBATCH --output=vg_minimap2_compare.%N.%j.out
#SBATCH --error=vg_minimap2_compare.%N.%j.err

########################################################################
# Carga de módulos
########################################################################
module load samtools
module load vg
module load minimap2

WORKDIR="/path/to/project/5.pangenome"
THREADS=20
mkdir -p "$WORKDIR"

########################################################################
# Listas de entrada (permite ejecución por subgrupos)
########################################################################
PART=${1:-1}

LIB_R1_LIST="/path/to/project/2.short_reads/sublists/R1_list_part_${PART}.txt"
LIB_R2_LIST="/path/to/project/2.short_reads/sublists/R2_list_part_${PART}.txt"
LIB_ID_LIST="/path/to/project/2.short_reads/sublists/ID_list_part_${PART}.txt"

echo "Ejecutando sublista: $PART"

mapfile -t R1_ARRAY < "$LIB_R1_LIST"
mapfile -t R2_ARRAY < "$LIB_R2_LIST"
mapfile -t ID_ARRAY < "$LIB_ID_LIST"

if [[ ${#R1_ARRAY[@]} -ne ${#R2_ARRAY[@]} || ${#R1_ARRAY[@]} -ne ${#ID_ARRAY[@]} ]]; then
    echo "Error: las listas de entrada no tienen el mismo número de elementos"
    exit 1
fi

########################################################################
# Preparación del grafo y referencia
########################################################################
COMBINED_GFA="$WORKDIR/merged_graph.gfa"
COMBINED_PREFIX="$WORKDIR/merged_graph"
COMBINED_GBZ="$COMBINED_PREFIX.gbz"
FIXED_GBZ="$COMBINED_PREFIX.giraffe.gbz"

COMBINED_REF="$WORKDIR/combined_reference.fa"
COMPARISON_FILE="$WORKDIR/comparison_summary.tsv"

if [[ ! -f "$COMBINED_GFA" ]]; then
    echo "Generando grafo combinado..."
    vg combine -p /path/to/gfa_files/*.gfa > "$COMBINED_GFA"
fi

if [[ ! -f "$COMBINED_GBZ" ]]; then
    echo "Indexando grafo para VG Giraffe..."
    vg autoindex -w giraffe --gfa "$COMBINED_GFA" --prefix "$COMBINED_PREFIX" --threads "$THREADS"

    if [[ ! -s "$FIXED_GBZ" ]]; then
        echo "Error: la indexación del grafo falló"
        exit 1
    fi

    mv "$FIXED_GBZ" "$COMBINED_GBZ"
fi

if [[ ! -f "$COMBINED_REF" ]]; then
    echo "Generando referencia lineal combinada..."
    cat $(cat /path/to/reference_list.txt) > "$COMBINED_REF"

    if [[ ! -s "$COMBINED_REF" ]]; then
        echo "Error: la referencia combinada está vacía"
        exit 1
    fi

    minimap2 -d "$WORKDIR/combined_reference.mmi" "$COMBINED_REF"
    samtools faidx "$COMBINED_REF"

    if [[ ! -f "$COMBINED_REF.fai" ]]; then
        echo "Error: no se pudo indexar la referencia"
        exit 1
    fi
fi

########################################################################
# Inicializar archivo de resultados
########################################################################
echo -e "Sample\tMethod\tMapped_Reads(%)\tAverage_MAPQ" > "$COMPARISON_FILE"

########################################################################
# Mapeo de lecturas
########################################################################
for i in "${!ID_ARRAY[@]}"; do

    LIB_ID="${ID_ARRAY[$i]}"
    LIB_R1="${R1_ARRAY[$i]}"
    LIB_R2="${R2_ARRAY[$i]}"

    PREFIX_VG="${LIB_ID}__vg"
    PREFIX_MINI="${LIB_ID}__minimap2"

    VG_BAM="$WORKDIR/${PREFIX_VG}.bam"
    MINIMAP2_BAM="$WORKDIR/${PREFIX_MINI}.bam"

    if [[ ! -f "$VG_BAM" ]]; then
        echo "Mapeo con VG (Giraffe): $LIB_ID"

        vg giraffe \
            -g "$COMBINED_GBZ" \
            -o gam \
            -p \
            -t "$THREADS" \
            -N "$LIB_ID" \
            -R "$LIB_ID" \
            -f "$LIB_R1" \
            -f "$LIB_R2" \
            > "$WORKDIR/${PREFIX_VG}.gam"

        vg surject \
            -x "$COMBINED_GBZ" \
            -b \
            -t "$THREADS" \
            "$WORKDIR/${PREFIX_VG}.gam" \
            > "$VG_BAM"

        rm -f "$WORKDIR/${PREFIX_VG}.gam"
    fi

    if [[ ! -f "$MINIMAP2_BAM" ]]; then
        echo "Mapeo con Minimap2: $LIB_ID"

        minimap2 \
            -t "$THREADS" \
            -ax sr \
            -R "@RG\tID:$LIB_ID\tSM:$LIB_ID" \
            "$COMBINED_REF" \
            "$LIB_R1" \
            "$LIB_R2" \
            | samtools view -bh > "$MINIMAP2_BAM"
    fi

done

echo "Proceso completado. Resultados en: $COMPARISON_FILE"
