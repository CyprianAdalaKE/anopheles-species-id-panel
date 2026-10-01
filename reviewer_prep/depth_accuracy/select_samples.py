#!/usr/bin/env python3
"""
Select ~8 samples per species spanning the native read-depth range
(quantile-stratified), for the controlled downsampling experiment.
"""
import csv
import numpy as np

EXCLUDE = {
    ("coluzzii", "barcode80"),
    ("coustani", "barcode84"),
    ("coustani", "barcode95"),
    ("arabiensis", "barcode06"),
}

MIN_READS = 5000  # need enough native reads to support the highest downsample tier

rows = []
with open("all_sample_readcounts.tsv") as f:
    for line in f:
        sp, bc, n, path = line.rstrip("\n").split("\t")
        if (sp, bc) in EXCLUDE:
            continue
        n = int(n)
        if n < MIN_READS:
            continue
        rows.append((sp, bc, n, path))

N_PER_SPECIES = 8
selected = []
for sp in sorted(set(r[0] for r in rows)):
    sub = sorted([r for r in rows if r[0] == sp], key=lambda r: r[2])
    if len(sub) <= N_PER_SPECIES:
        selected.extend(sub)
        continue
    idxs = np.linspace(0, len(sub) - 1, N_PER_SPECIES).astype(int)
    idxs = sorted(set(idxs))
    for i in idxs:
        selected.append(sub[i])

with open("selected_samples.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["species", "barcode", "native_reads", "fastq_path"])
    w.writerows(selected)

print(f"Selected {len(selected)} samples")
for sp, bc, n, path in selected:
    print(f"  {sp:12s} {bc:12s} n={n}")
