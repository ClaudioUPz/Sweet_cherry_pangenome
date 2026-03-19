# Step 5: Core, shared, and private genome analysis

This step describes how core, shared, and private genomic regions were quantified from the pangenome graph and visualized per genome/haplotype.

---

## Overview

The workflow includes:

1. Parsing the pangenome graph (GFA)
2. Assigning graph nodes to genomes using path information
3. Classifying sequences as core, shared, or private
4. Generating per-genome statistics
5. Visualizing results

---

## Requirements

- Python 3
- R
- Required Python libraries:
  - argparse
  - gzip
- Required R packages:
  - readr
  - dplyr
  - tidyr
  - ggplot2
  - scales

---

## Step 1: Compute core, shared, and private sequences

A custom Python script was used to classify graph nodes based on their presence across genomes.

### Description

- Each node in the GFA graph is assigned to:
  - **Core**: present in ≥ 90% of genomes (default threshold)
  - **Shared**: present in more than one but not core
  - **Private**: present in only one genome

- The output is a per-sample summary of base pairs (bp) in each category.

### Script

Save as `compute_core_shared_private.py`:

```python
#!/usr/bin/env python3
# (script content unchanged, truncated here for readability in README)

Usage
python compute_core_shared_private.py \
  merged_ReginaC.full.gfa \
  tree_path.tsv \
  per_sample_core_shared_private.tsv \
  --core_fraction 0.9
Input

merged_ReginaC.full.gfa – pangenome graph in GFA format

tree_path.tsv – mapping of paths to genome/sample names

Output

per_sample_core_shared_private.tsv with columns:

sample    total_bp    core_bp    shared_bp    private_bp
Step 2: Plot core, shared, and private fractions

An R script was used to generate bar plots showing the fraction of genomic bases per category.

Script

Save as plot_core_shared_private.R:

#!/usr/bin/env Rscript
# (script content unchanged, truncated here for readability in README)

Usage
Rscript plot_core_shared_private.R
Output

Fig1C_core_shared_private.pdf

Fig1C_core_shared_private.png

Notes

The core_fraction parameter defines the minimum proportion of genomes required for a node to be considered core (default = 0.9).

Sample names are standardized before plotting to ensure consistent labeling.

This analysis is based on graph topology and path assignments rather than linear genome coordinates.

Recommended directory structure
scripts/05_core_shared_private/
├── README.md
├── compute_core_shared_private.py
└── plot_core_shared_private.R
