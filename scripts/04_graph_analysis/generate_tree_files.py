#!/usr/bin/env python3
import re
from ete3 import Tree

path_file = "path_ReginaC.txt"
tree_file = "genome_kmer_distance.tree"

out_order = "tree_order.txt"
out_path  = "tree_path.tsv"

def normalize(name: str) -> str:
    name = re.sub(r"hap1", "h1", name)
    name = re.sub(r"hap2", "h2", name)
    return name

def chr_sort_key(p: str) -> int:
    m = re.search(r"chr(\d+)", p)
    return int(m.group(1)) if m else 999

t = Tree(tree_file)
genomes_ordered = [normalize(leaf.name) for leaf in t.get_leaves()]

with open(path_file) as f:
    paths = [x.strip() for x in f if x.strip()]

paths_by_genome = {}
for p in paths:
    base = re.split(r"#0#|__", p)[0]
    paths_by_genome.setdefault(base, []).append(p)

with open(out_order, "w") as out:
    for genome in genomes_ordered:
        if genome not in paths_by_genome:
            print(f"[WARN] Genome {genome} not found in paths")
            continue
        for p in sorted(paths_by_genome[genome], key=chr_sort_key):
            out.write(p + "\n")

with open(out_path, "w") as out:
    out.write("path_name\tsample_name\n")
    for genome in genomes_ordered:
        if genome not in paths_by_genome:
            continue
        for p in sorted(paths_by_genome[genome], key=chr_sort_key):
            out.write(f"{p}\t{genome}\n")

