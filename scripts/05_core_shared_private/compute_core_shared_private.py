#!/usr/bin/env python3
import sys
import math
import argparse
import gzip
from collections import defaultdict

def open_maybe_gzip(path, mode="rt"):
    """
    Abre un archivo normal o .gz de forma transparente.
    mode por defecto es texto ("rt").
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode)
    else:
        return open(path, mode=mode, encoding="utf-8")

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calcular, para cada sample/genoma, cuánto BP pertenece a regiones "
            "core / shared / private usando un GFA de pangenoma y un tree_path TSV."
        )
    )
    parser.add_argument("gfa", help="GFA del pangenoma (por ejemplo merged_Regina.full.gfa)")
    parser.add_argument("tree_path_tsv", help="TSV con columnas path_name y sample_name")
    parser.add_argument("out_tsv", help="TSV de salida con resumen por sample")
    parser.add_argument(
        "--core_fraction",
        type=float,
        default=0.9,
        help="Fracción mínima de samples en la que debe estar un nodo para ser CORE (default 0.9)"
    )
    return parser.parse_args()

def load_tree_path(tree_path_file):
    """
    Lee el TSV (tree_path_full.fixed.tsv) y devuelve:
      - dict_path_to_sample: path_name -> sample_name
      - sample_list: lista de samples en el orden en que aparecen
    """
    dict_path_to_sample = {}
    sample_list = []
    seen_samples = set()

    with open(tree_path_file, "r", encoding="utf-8") as f:
        header = f.readline()
        header = header.strip().split("\t")
        if len(header) < 2 or header[0] != "path_name":
            sys.stderr.write(
                "Advertencia: primera línea de tree_path no parece header estándar "
                "(path_name\tsample_name). Se usará igualmente.\n"
            )

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            path_name, sample_name = parts[0], parts[1]
            dict_path_to_sample[path_name] = sample_name
            if sample_name not in seen_samples:
                seen_samples.add(sample_name)
                sample_list.append(sample_name)

    return dict_path_to_sample, sample_list

def canonicalize_path_name(path_name_raw, dict_path_to_sample):
    """
    Intenta mapear el nombre del path tal como sale del GFA
    a uno presente en el tree_path TSV, tolerando un '#0' final.
    """
    # Intento 1: tal cual
    if path_name_raw in dict_path_to_sample:
        return path_name_raw

    # Intento 2: si termina en "#0", probamos sin ese sufijo
    if path_name_raw.endswith("#0"):
        alt = path_name_raw[:-2]
        if alt in dict_path_to_sample:
            return alt

    # Si no encontramos correspondencia, devolvemos None
    return None

def main():
    args = parse_args()

    gfa_path = args.gfa
    tree_path_tsv = args.tree_path_tsv
    out_tsv = args.out_tsv
    core_fraction = args.core_fraction

    sys.stderr.write(f"Usando GFA: {gfa_path}\n")
    sys.stderr.write(f"Usando tree_path: {tree_path_tsv}\n")
    sys.stderr.write(f"core_fraction = {core_fraction}\n")

    # 1) Cargar tree_path TSV
    dict_path_to_sample, sample_list = load_tree_path(tree_path_tsv)

    n_samples = len(sample_list)
    sys.stderr.write(f"Número de samples (genomas/haplotipos): {n_samples}\n")
    for i, s in enumerate(sample_list):
        sys.stderr.write(f"  [{i}] {s}\n")

    # Mapear sample -> índice
    sample_to_idx = {s: i for i, s in enumerate(sample_list)}

    # 2) Recorrer GFA: obtener longitud de cada nodo (S-lines)
    #    y, al mismo tiempo, construir la máscara de presencia por nodo (P-lines)
    node_len = {}
    node_mask = defaultdict(int)  # node_id -> bitmask de samples
    paths_in_gfa = set()
    used_paths = set()

    total_nodes = 0
    total_paths = 0

    with open_maybe_gzip(gfa_path, "rt") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            line = line.rstrip("\n")
            if not line:
                continue

            if line[0] == "S":
                # Formato S: S <node_id> <seq> ...
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                node_id = parts[1]
                seq = parts[2]
                node_len[node_id] = len(seq)
                total_nodes += 1

            elif line[0] == "P":
                # Formato P: P <path_name> <segments> <overlaps> ...
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                path_name_raw = parts[1]
                segments = parts[2]

                paths_in_gfa.add(path_name_raw)
                total_paths += 1

                canonical = canonicalize_path_name(path_name_raw, dict_path_to_sample)
                if canonical is None:
                    # path no está en el tree_path, lo ignoramos
                    continue

                sample_name = dict_path_to_sample[canonical]
                if sample_name not in sample_to_idx:
                    # Esto no debería pasar, pero por si acaso
                    continue
                s_idx = sample_to_idx[sample_name]
                used_paths.add(canonical)

                # segments: lista tipo "12+,13-,14+"
                for tok in segments.split(","):
                    tok = tok.strip()
                    if not tok:
                        continue
                    # Nodo es todo menos el último carácter (+ o -)
                    node_id = tok[:-1]
                    # Marcamos el bit del sample
                    node_mask[node_id] |= (1 << s_idx)

    sys.stderr.write(f"Nodos leídos (S-lines): {total_nodes}\n")
    sys.stderr.write(f"Paths en GFA (P-lines): {total_paths}\n")
    sys.stderr.write(f"Paths usados (presentes en tree_path): {len(used_paths)}\n")

    # Contar cuántos nodos tienen al menos un sample
    nodes_with_mask = sum(1 for m in node_mask.values() if m != 0)
    sys.stderr.write(f"Nodos con máscara de presencia: {nodes_with_mask}\n")

    # 3) Calcular umbral de CORE en número de samples
    if 0 < core_fraction <= 1.0:
        core_k = math.ceil(core_fraction * n_samples)
    else:
        # Si el usuario pone un número >= 1, lo interpretamos como número de samples
        core_k = int(core_fraction)
    sys.stderr.write(f"Umbral CORE: nodo presente en >= {core_k} samples\n\n")

    # 4) Inicializar contadores por sample
    total_bp = [0] * n_samples
    core_bp = [0] * n_samples
    shared_bp = [0] * n_samples
    private_bp = [0] * n_samples

    # Función helper para contar bits
    def popcount(x):
        return x.bit_count()  # Python 3.8+

    # 5) Recorrer nodos y clasificar
    processed = 0
    for node_id, mask in node_mask.items():
        if mask == 0:
            continue
        length = node_len.get(node_id, 0)
        if length == 0:
            continue

        k = popcount(mask)

        # Clasificación del nodo según cuántos samples lo tienen
        if k >= core_k:
            category = "core"
        elif k == 1:
            category = "private"
        else:
            category = "shared"

        # Recorremos los bits activos del mask
        m = mask
        while m:
            lowest_bit = m & -m
            s_idx = (lowest_bit.bit_length() - 1)
            total_bp[s_idx] += length
            if category == "core":
                core_bp[s_idx] += length
            elif category == "shared":
                shared_bp[s_idx] += length
            else:  # private
                private_bp[s_idx] += length
            m ^= lowest_bit

        processed += 1
        # Puedes descomentar esto si quieres ver progreso (pero puede ser mucho output)
        # if processed % 1000000 == 0:
        #     sys.stderr.write(f"  Procesados {processed} nodos con máscara\n")

    # 6) Escribir salida
    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("sample\ttotal_bp\tcore_bp\tshared_bp\tprivate_bp\n")
        for i, sample in enumerate(sample_list):
            out.write(
                f"{sample}\t{total_bp[i]}\t{core_bp[i]}"
                f"\t{shared_bp[i]}\t{private_bp[i]}\n"
            )

    sys.stderr.write(f"Escribí resumen por sample en: {out_tsv}\n\n")
    sys.stderr.write("Pequeño resumen:\n")
    for i, sample in enumerate(sample_list):
        sys.stderr.write(
            f"{sample}: total={total_bp[i]} bp | "
            f"Core={core_bp[i]} | Shared={shared_bp[i]} | Private={private_bp[i]}\n"
        )

if __name__ == "__main__":
    main()
