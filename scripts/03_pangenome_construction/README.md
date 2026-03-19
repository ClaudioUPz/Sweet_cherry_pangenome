# Step 3: Pangenome construction with MiniCactus

This step describes how genome assemblies were prepared and used to construct a chromosome-level pangenome with MiniCactus.

## Overview

The workflow includes:

1. Splitting genome assemblies by chromosome
2. Preparing chromosome-specific `.dataset` files
3. Running MiniCactus with Singularity
4. Validating output graphs with `vg`

---

## Requirements

- Python 3
- Biopython
- Singularity
- MiniCactus
- `vg`
- SLURM-based HPC environment

---

## Step 1: Split genome assemblies by chromosome

Genome FASTA files were split into separate chromosome-specific FASTA files. Each chromosome sequence was written to its corresponding chromosome directory.

### Script

Save the following script as `split_by_chr.py`:

```python
import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import csv

def split_genomes_by_chromosome_and_log():
    base_dest = "/home/curra/pangenome/5.pangenome"
    chromosomes = [f"chr{i}" for i in range(1, 9)]
    fasta_files = sorted([f for f in os.listdir(".") if f.endswith(".fasta")])

    for chr_id in chromosomes:
        os.makedirs(os.path.join(base_dest, chr_id), exist_ok=True)

    log_records = []

    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            id_split = record.id.split("__")
            if len(id_split) == 2 and id_split[0] in chromosomes:
                chr_id = id_split[0]
                genome = id_split[1]
                output_path = os.path.join(base_dest, chr_id, f"{genome}_{chr_id}.fa")

                with open(output_path, "w") as out:
                    writer = FastaWriter(out, wrap=None)
                    writer.write_header()
                    writer.write_record(record)
                    writer.write_footer()

                log_records.append({
                    "source_file": fasta_file,
                    "chromosome": chr_id,
                    "genome_name": genome,
                    "output_file": output_path
                })
                print(f"Saved: {output_path}")

    log_path = os.path.join(base_dest, "fasta_split_log.tsv")
    with open(log_path, "w", newline='') as log_file:
        writer = csv.DictWriter(
            log_file,
            delimiter='\t',
            fieldnames=["source_file", "chromosome", "genome_name", "output_file"]
        )
        writer.writeheader()
        for row in log_records:
            writer.writerow(row)

    print(f"Log written to: {log_path}")

if __name__ == "__main__":
    split_genomes_by_chromosome_and_log()
Usage

Run the script in a directory containing all preprocessed genome FASTA files:

python split_by_chr.py
Output

The script creates chromosome-specific FASTA files in directories such as:

/home/curra/pangenome/5.pangenome/chr1/
/home/curra/pangenome/5.pangenome/chr2/
...
/home/curra/pangenome/5.pangenome/chr8/

It also generates a log file:

/home/curra/pangenome/5.pangenome/fasta_split_log.tsv
Step 2: Prepare .dataset files for MiniCactus

MiniCactus requires a chromosome-specific .dataset file.

Structure of the .dataset file

Each .dataset file contains:

A first line with the phylogenetic tree in Newick format

One line per genome, with:

the genome name used in the tree

the full path to the corresponding chromosome FASTA file

The genome entries must be listed in the same order as they appear in the Newick tree.

Example
(ReginaC:0.00225474,(ReginaF:0.0059228,(V2775_h2:0.00340396,...):0.00011867);

ReginaC /home/curra/pangenome/5.pangenome/chr2/ReginaC_chr2.fa
ReginaF /home/curra/pangenome/5.pangenome/chr2/ReginaF_chr2.fa
V2775_h2 /home/curra/pangenome/5.pangenome/chr2/V2775_h2_chr2.fa
V3382_h1 /home/curra/pangenome/5.pangenome/chr2/V3382_h1_chr2.fa
V2076_h1 /home/curra/pangenome/5.pangenome/chr2/V2076_h1_chr2.fa
...
Important notes

Genome names in the Newick tree must match the names used in the dataset file.

Genome names should not contain dots (.); use underscores (_) instead.

A separate .dataset file must be created for each chromosome.

Consistent naming between tree labels, FASTA filenames, and FASTA headers is required.

Reference naming note

For compatibility with downstream analyses, the reference haplotype name was simplified from Regina_h1 to ReginaC, and the second haplotype was renamed to ReginaF.

Step 3: Build the MiniCactus container and run the pangenome workflow
Build the Singularity image
singularity build cactus.sif docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13
Export the container path
export CACTUS_SIF=/home/curra/cactus.sif
SLURM job script

Save the following script as run_minicactus.sh:

#!/bin/bash
#SBATCH --job-name=cactus_run
#SBATCH --time=96:00:00
#SBATCH -o cactus_chr%a_sakuromics.%j.out
#SBATCH -e cactus_chr%a_sakuromics.%j.err
#SBATCH --mem=512G
#SBATCH --cpus-per-task=48
#SBATCH --array=1-8

module purge

CHR_ID=${SLURM_ARRAY_TASK_ID}
WORKDIR="/home/curra/pangenome/5.pangenome"
DATASET="${WORKDIR}/Pavium_chr${CHR_ID}.dataset"
ANALYSIS_ID="chr${CHR_ID}_run"
REFERENCE="ReginaC"
VCF_REFS="ReginaC"

DIR_WORKDIR="${WORKDIR}/chr${CHR_ID}"
DIR_CACTUS_OUTPUT="${DIR_WORKDIR}/output"
DIR_CACTUS_SCRATCH="${DIR_WORKDIR}/scratch"
DIR_CACTUS_JS="${DIR_WORKDIR}/js"

mkdir -p "$DIR_CACTUS_OUTPUT" "$DIR_CACTUS_SCRATCH"
rm -rf "$DIR_CACTUS_JS"

singularity exec -B ${PWD}:${PWD} /home/curra/cactus.sif \
    cactus-pangenome "$DIR_CACTUS_JS" "$DATASET" \
    --outDir "$DIR_CACTUS_OUTPUT" \
    --logFile "$DIR_CACTUS_OUTPUT/run.log" \
    --workDir "$DIR_CACTUS_SCRATCH" \
    --outName "$ANALYSIS_ID" \
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
Usage

Submit the SLURM array job:

sbatch run_minicactus.sh
Output

For each chromosome, MiniCactus generates graph-based pangenome outputs in:

/home/curra/pangenome/5.pangenome/chrX/output/
Step 4: Validate MiniCactus graphs with vg

The .d1.vg graph files were validated using vg v1.65.0.

SLURM job script

Save the following script as validate_vg_graphs.sh:

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

BASE_DIR="/home/curra/pangenome/5.pangenome"

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
Usage
sbatch validate_vg_graphs.sh
