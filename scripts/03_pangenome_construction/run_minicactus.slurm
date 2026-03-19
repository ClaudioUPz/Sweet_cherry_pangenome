#!/bin/bash
#SBATCH --job-name=cactus_run
#SBATCH --time=96:00:00
#SBATCH -o cactus_chr%a_sakuromics.%j.out
#SBATCH -e cactus_chr%a_sakuromics.%j.err
#SBATCH --mem=512G
#SBATCH --cpus-per-task=48
#SBATCH --array=1-8

# Configuración del entorno
module purge

# Definir variables
CHR_ID=${SLURM_ARRAY_TASK_ID}
WORKDIR="/path/to/workdir"
DATASET="${WORKDIR}/Pavium_chr${CHR_ID}.dataset"
ANALYSE_ID="chr${CHR_ID}_run"
REFERENCE="Regina"
VCF_REFS="Regina"

# Directorios de salida
DIR_WORKDIR="${WORKDIR}/chr${CHR_ID}"
DIR_CACTUS_OUTPUT="${DIR_WORKDIR}/output"
DIR_CACTUS_SCRATCH="${DIR_WORKDIR}/scratch"
DIR_CACTUS_JS="${DIR_WORKDIR}/js"

# Crear carpetas necesarias (o limpiar si ya existen)
mkdir -p "$DIR_CACTUS_OUTPUT" "$DIR_CACTUS_SCRATCH"
rm -rf "$DIR_CACTUS_JS"

# Ejecutar MiniCactus con VCF habilitado
singularity exec -B ${PWD}:${PWD} /home/curra/cactus.sif \
    cactus-pangenome "$DIR_CACTUS_JS" "$DATASET" \
    --outDir "$DIR_CACTUS_OUTPUT" \
    --logFile "$DIR_CACTUS_OUTPUT/run.log" \
    --workDir "$DIR_CACTUS_SCRATCH" \
    --outName "$ANALYSE_ID" \
    --reference "$REFERENCE" \
    --vcf --vcfReference "$VCF_REFS" \
    --filter 1 \
    --haplo \
    --viz \
    --odgi \
    --giraffe clip filter full \
    --chrom-vg clip filter full \
    --chrom-og \
    --gbz clip filter full \
    --gfa clip filter full \
    --consCores 48 \
    --mgMemory 450Gi \
    --indexMemory 450Gi \
    --indexCores 48
