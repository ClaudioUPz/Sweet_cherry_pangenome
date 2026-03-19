# Step 1: Sequence preprocessing

This step includes preprocessing of genome FASTA files before downstream comparative analyses. It covers:

1. FASTA header renaming
2. FASTA sequence reformatting to single-line entries

These formatting steps help standardize genome assemblies for subsequent analyses.

---

## 1. Rename genome FASTA headers

This script standardizes FASTA headers across genome assemblies to ensure consistent naming.

### Description

The script renames chromosome headers in FASTA files using the following format:

```text
chrX__GenomeName_hY

Where:

chrX = original chromosome name

GenomeName = derived from the FASTA filename

hY = haplotype identifier (h1 or h2, if present)

Dots (.) in filenames are replaced with underscores (_) to avoid formatting issues.

Requirements

Python 3

Biopython

Install Biopython if needed:

pip install biopython
Script

Save the following script as rename_headers.py:

import os
from Bio import SeqIO

def get_genome_name_and_haplotype(filename):
    base = filename.replace(".fasta", "")
    base = base.replace(".", "_")
    if "_hap1" in base:
        return base.replace("_hap1", ""), "h1"
    elif "_hap2" in base:
        return base.replace("_hap2", ""), "h2"
    else:
        return base, None

def rename_headers(input_fasta):
    genome, hap = get_genome_name_and_haplotype(os.path.basename(input_fasta))
    output_fasta = input_fasta.replace(".fasta", "_renamed.fasta")

    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            base_id = record.id.split()[0]
            if base_id.startswith("chr"):
                new_id = f"{base_id}__{genome}"
                if hap:
                    new_id += f"_{hap}"
                record.id = new_id
                record.description = ""
            SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    for filename in os.listdir("."):
        if filename.endswith(".fasta"):
            print(f"Processing: {filename}")
            rename_headers(filename)
Usage

Place all FASTA files in the same directory as the script and run:

python rename_headers.py
Input

FASTA files with chromosome headers (e.g., chr1, chr2, etc.)

Output

Renamed FASTA files with suffix:

*_renamed.fasta
Notes

Only headers starting with chr are modified.

Haplotype information is automatically extracted from filenames containing _hap1 or _hap2.

2. Reformat FASTA sequences to single-line format

This script rewrites FASTA files so that each sequence is written on a single line instead of multiple fixed-width lines.

Description

Some downstream tools require or work better with FASTA files where each sequence is stored on a single line. This script reformats all .fasta files in the working directory accordingly.

Requirements

Python 3

Biopython

Script

Save the following script as reformat_fasta_single_line.py:

import os
from Bio import SeqIO

def rewrite_fasta_single_line(input_fasta):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    with open(input_fasta, "w") as out:
        for record in records:
            out.write(f">{record.id}\n{str(record.seq)}\n")
    print(f"Reformatted file: {input_fasta}")

def process_all_fasta():
    for filename in os.listdir("."):
        if filename.endswith(".fasta"):
            rewrite_fasta_single_line(filename)

if __name__ == "__main__":
    process_all_fasta()
Usage

Place all FASTA files in the same directory as the script and run:

python reformat_fasta_single_line.py
Input

FASTA files in standard multi-line format

Output

The original FASTA files are overwritten with single-line sequence format

Notes

This script modifies files in place.

It is recommended to keep a backup copy of the original FASTA files before running this step.
