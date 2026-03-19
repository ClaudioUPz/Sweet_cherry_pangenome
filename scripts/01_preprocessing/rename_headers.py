import os
import re
from Bio import SeqIO

def obtener_nombre_genoma_y_haplotipo(nombre_archivo):
    # Elimina la extensión .fasta
    base = nombre_archivo.replace(".fasta", "")
    # Reemplaza "." por "_" (para evitar problemas)
    base = base.replace(".", "_")
    # Extrae haplotipo si está en el nombre
    if "_hap1" in base:
        return base.replace("_hap1", ""), "h1"
    elif "_hap2" in base:
        return base.replace("_hap2", ""), "h2"
    else:
        return base, None

def renombrar_headers(input_fasta):
    genoma, hap = obtener_nombre_genoma_y_haplotipo(os.path.basename(input_fasta))
    output_fasta = input_fasta.replace(".fasta", "_renombrado.fasta")

    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            header_parts = record.id.split()
            base_id = header_parts[0]
            if base_id.startswith("chr"):
                nuevo_id = f"{base_id}__{genoma}"
                if hap:
                    nuevo_id += f"_{hap}"
                record.id = nuevo_id
                record.description = ""
            SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    for archivo in os.listdir("."):
        if archivo.endswith(".fasta"):
            print(f"Procesando: {archivo}")
            renombrar_headers(archivo)
