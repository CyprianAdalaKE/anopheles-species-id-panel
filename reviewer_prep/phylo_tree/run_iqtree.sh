#!/bin/bash
# Maximum-likelihood tree from the SNP alignment written by build_alignment.py.
# Command recovered from iqtree_run.log of the original run (~73 min wall-clock).
set -euo pipefail
cd "$(dirname "$0")"
iqtree -s snp_alignment.fasta -st DNA -m GTR+G -B 1000 -T AUTO -ntmax 16 -pre snp_tree -redo
