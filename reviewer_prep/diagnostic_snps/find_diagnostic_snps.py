#!/usr/bin/env python3
"""
Identify the strongest candidate diagnostic SNP per locus: for every SNP,
compute allele frequency in each of the 6 species and find the species pair
with the largest frequency difference. Report the single best SNP per locus
(by max |delta freq| across any pair) -- a practical "if you had to pick one
site per locus to discriminate two species, this is it" table.
"""
import csv
import itertools
import numpy as np

PCA_DIR = ".."  # run from anopheles-panel/pca_validation
BED = "../nomadsSID.amplicons.bed"

SPECIES_LIST = ["gambiae", "coluzzii", "arabiensis", "funestus", "stephensi", "coustani"]


def gt_to_dosage(gt):
    if gt in (".", "./.", ".|."):
        return np.nan
    sep = "/" if "/" in gt else "|"
    a, b = gt.split(sep)
    if a == "." or b == ".":
        return np.nan
    return int(a) + int(b)


with open("manifest.csv") as f:
    manifest = {row["sample_id"]: row["species"] for row in csv.DictReader(f)}
with open("results/sample_order.txt") as f:
    sample_order = [l.strip() for l in f if l.strip()]


def path_to_sample_id(p):
    parts = p.split("/")
    return f"{parts[-2]}_{parts[-1].split('.')[0]}"


sample_ids = np.array([path_to_sample_id(s) for s in sample_order])
species = np.array([manifest[s] for s in sample_ids])

# drop the known high-missingness coluzzii batch, same QC as before
sites, chroms, poss, refs, alts = [], [], [], [], []
rows = []
with open("results/genotype_matrix.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = parts[0:4]
        chroms.append(chrom); poss.append(int(pos)); refs.append(ref); alts.append(alt)
        rows.append([gt_to_dosage(g) for g in parts[4:]])
G = np.array(rows, dtype=float).T  # samples x SNPs
chroms = np.array(chroms); poss = np.array(poss); refs = np.array(refs); alts = np.array(alts)

miss = np.isnan(G).mean(axis=1)
keep_samples = miss <= 0.5
G = G[keep_samples]
species_k = species[keep_samples]

loci = []
with open(BED) as f:
    for line in f:
        c, s, e, name = line.strip().split("\t")
        loci.append((name, c, int(s), int(e)))

pairs = list(itertools.combinations(SPECIES_LIST, 2))

results = []
for name, c, s, e in loci:
    site_idx = np.where((chroms == c) & (poss >= s) & (poss <= e))[0]
    best = None
    for j in site_idx:
        col = G[:, j]
        freqs = {}
        for sp in SPECIES_LIST:
            m = species_k == sp
            vals = col[m]
            vals = vals[~np.isnan(vals)]
            if len(vals) < 5:
                freqs[sp] = np.nan
            else:
                freqs[sp] = vals.mean() / 2.0
        for spA, spB in pairs:
            if np.isnan(freqs[spA]) or np.isnan(freqs[spB]):
                continue
            d = abs(freqs[spA] - freqs[spB])
            if best is None or d > best["delta"]:
                best = {
                    "locus": name, "chrom": c, "pos": poss[j], "ref": refs[j], "alt": alts[j],
                    "spA": spA, "spB": spB, "freqA": freqs[spA], "freqB": freqs[spB], "delta": d,
                    "freqs_all": dict(freqs),
                }
    if best:
        results.append(best)

results.sort(key=lambda r: -r["delta"])
print(f"{'Locus':12s} {'Pos':>10s} {'Ref>Alt':8s} {'Pair':22s} {'freqA':>6s} {'freqB':>6s} {'|delta|':>7s}")
for r in results:
    print(f"{r['locus']:12s} {r['pos']:10d} {r['ref']}>{r['alt']:6s} {r['spA']}-{r['spB']:14s} "
          f"{r['freqA']:6.2f} {r['freqB']:6.2f} {r['delta']:7.2f}")

with open("../reviewer_prep/diagnostic_snps/diagnostic_snps.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Locus", "Position (AgamP4)", "Ref>Alt", "Best-discriminated pair",
                "Freq in sp. A", "Freq in sp. B", "|Freq diff|"] + [f"freq_{sp}" for sp in SPECIES_LIST])
    for r in results:
        w.writerow([r["locus"], r["pos"], f"{r['ref']}>{r['alt']}", f"{r['spA']} vs {r['spB']}",
                    f"{r['freqA']:.2f}", f"{r['freqB']:.2f}", f"{r['delta']:.2f}"] +
                   [f"{r['freqs_all'][sp]:.2f}" if not np.isnan(r['freqs_all'][sp]) else "NA" for sp in SPECIES_LIST])
print("\nSaved diagnostic_snps.tsv")
