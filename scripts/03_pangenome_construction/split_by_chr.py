import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import csv

def split_genomes_by_chromosome_y_log():
    out_dir = "/path/to/outdir"
    chromosomes = [f"chr{i}" for i in range(1, 9)]
    files = sorted([f for f in os.listdir(".") if f.endswith(".fasta")])

    # Make directories chr1/ ... chr8/ if they don't exist
    for chr_id in chromosomes:
        os.makedirs(os.path.join(out_dir, chr_id), exist_ok=True)

    log_register = []

    for file in files:
        for record in SeqIO.parse(file, "fasta"):
            id_split = record.id.split("__")
            if len(id_split) == 2 and id_split[0] in chromosomes:
                chr_id = id_split[0]
                genome = id_split[1]
                output_path = os.path.join(out_dir, chr_id, f"{genome}_{chr_id}.fa")
                
                with open(output_path, "w") as out:
                    writer = FastaWriter(out, wrap=None)  # wrap=None → one line
                    writer.write_header()
                    writer.write_record(record)
                    writer.write_footer()

                log_register.append({
                    "original_file": file,
                    "chromosome": chr_id,
                    "genome_name": genome,
                    "output_file": output_path
                })
                print(f"✔ Saved: {output_path}")

    # Write log
    log_path = os.path.join(out_dir, "fasta_division_register.tsv")
    with open(log_path, "w", newline='') as log_file:
        writer = csv.DictWriter(log_file, delimiter='\t', fieldnames=["original_file", "chromosome", "genome_name", "output_file"])
        writer.writeheader()
        for row in log_register:
            writer.writerow(row)

    print(f"✔ Log generated: {log_path}")

if __name__ == "__main__":
    split genomes by chromosome_y_log()
