# Step 4: Graph analysis and pangenome statistics

This step describes how graph-based statistics were computed from the pangenome and how pangenome growth was estimated using Panacus.

---

## Overview

The workflow includes:

1. Indexing graph files
2. Extracting paths and regions
3. Computing graph statistics
4. Merging chromosome graphs
5. Preparing inputs for Panacus
6. Estimating pangenome growth

---

## Requirements

- vg (v1.65.0)
- Python 3
- ete3 (Python library)
- panacus
- sed

---

## Step 1: Index graph files

Graph files (`.vg`) were indexed using `vg index`:

```bash
vg index -x chr1__ReginaC.xg chr1__ReginaC.vg
Step 2: Extract graph paths

Paths were extracted from the index:

vg paths -x chr1__ReginaC.xg -L

This allows identification of path names, including the reference genome path.

Step 3: Compute graph statistics

Basic statistics were computed for each chromosome graph:

vg stats chr1_run.full.gbz -z -N -E -l > chr1_run.full_stats.txt
Step 4: Merge chromosome graphs

Chromosome-level graphs were combined into a single genome-wide graph:

vg combine \
chr1/backup/output/chr1_run.chroms/chr1__ReginaC.full.vg \
chr2/backup/output/chr2_run.chroms/chr2__ReginaC.full.vg \
chr3/backup/output/chr3_run.chroms/chr3__ReginaC.full.vg \
chr4/backup/output/chr4_run.chroms/chr4__ReginaC.full.vg \
chr5/backup/output/chr5_run.chroms/chr5__ReginaC.full.vg \
chr6/backup/output/chr6_run.chroms/chr6__ReginaC.full.vg \
chr7/backup/output/chr7_run.chroms/chr7__ReginaC.full.vg \
chr8/backup/output/chr8_run.chroms/chr8__ReginaC.full.vg \
> merged_ReginaC.full.vg
Step 5: Convert graph to GFA format
vg view -g merged_ReginaC.full.vg > merged_ReginaC.full.gfa
Step 6: Extract paths from merged graph
vg paths -x merged_ReginaC.full.vg -L

The resulting file contains all paths in the graph and is used for downstream analysis.

Step 7: Prepare input files for Panacus

A custom script was used to generate the required files:

tree_order.txt

tree_path.tsv

These files define:

the order of genomes

the mapping between graph paths and genomes

Script

Save as generate_tree_files.py:

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
Post-processing

Remove trailing #0 from path names:

sed -E 's/__([A-Za-z0-9_]+)#0$/__\1/' tree_order.txt > tree_order.fixed.txt
Step 8: Compute pangenome growth with Panacus

Activate the Panacus environment:

conda activate panacus

Run the ordered histogram growth analysis:

panacus ordered-histgrowth merged_ReginaC.full.gfa \
  --quorum 1,0 \
  --count bp \
  --groupby tree_path.tsv \
  --order tree_order.fixed.txt \
  -t 8 > merged_growth_size_ReginaC.tsv
Output

merged_ReginaC.full.vg – merged graph

merged_ReginaC.full.gfa – graph in GFA format

tree_order.txt – genome order

tree_path.tsv – path-to-genome mapping

merged_growth_size_ReginaC.tsv – pangenome growth statistics
