import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import csv

def dividir_genomas_por_cromosoma_y_log():
    base_destino = "/home/curra/pangenome/5.pangenome"
    cromosomas = [f"chr{i}" for i in range(1, 9)]
    archivos = sorted([f for f in os.listdir(".") if f.endswith(".fasta")])

    # Crear carpetas chr1/ ... chr8/ si no existen
    for chr_id in cromosomas:
        os.makedirs(os.path.join(base_destino, chr_id), exist_ok=True)

    log_registros = []

    for archivo in archivos:
        for record in SeqIO.parse(archivo, "fasta"):
            id_split = record.id.split("__")
            if len(id_split) == 2 and id_split[0] in cromosomas:
                chr_id = id_split[0]
                genoma = id_split[1]
                output_path = os.path.join(base_destino, chr_id, f"{genoma}_{chr_id}.fa")
                
                with open(output_path, "w") as out:
                    writer = FastaWriter(out, wrap=None)  # wrap=None → una sola línea
                    writer.write_header()
                    writer.write_record(record)
                    writer.write_footer()

                log_registros.append({
                    "archivo_origen": archivo,
                    "cromosoma": chr_id,
                    "nombre_genoma": genoma,
                    "archivo_salida": output_path
                })
                print(f"✔ Guardado: {output_path}")

    # Escribir log
    log_path = os.path.join(base_destino, "registro_division_fasta.tsv")
    with open(log_path, "w", newline='') as log_file:
        writer = csv.DictWriter(log_file, delimiter='\t', fieldnames=["archivo_origen", "cromosoma", "nombre_genoma", "archivo_salida"])
        writer.writeheader()
        for row in log_registros:
            writer.writerow(row)

    print(f"✔ Log generado: {log_path}")

if __name__ == "__main__":
    dividir_genomas_por_cromosoma_y_log()
