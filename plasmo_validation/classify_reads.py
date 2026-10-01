#!/usr/bin/env python3
"""
Competitive alignment-based Plasmodium/host classification.

For each infectivity sample, reads were aligned separately to (a) the AgamP4
host genome and (b) a combined 4-species Plasmodium reference (falciparum,
vivax, malariae, ovale). This script compares each read's best alignment on
each side (percent identity x query coverage) and classifies it as
Plasmodium-positive, host-positive, or unclassified/ambiguous -- replacing
Kraken2/Bracken's k-mer/LCA call with a direct, whole-read alignment
comparison, which should be far less susceptible to the shared-conserved-
k-mer problem identified for the 18S target.
"""
import csv
import glob
import os
import pysam

BAM_DIR = "bams"
RESULTS_DIR = "results"

IDENTITY_THRESH = 0.85
COVERAGE_THRESH = 0.70
MARGIN = 0.03  # winning side's (identity*coverage) score must exceed the other by this much

GROUPS = {
    **{f"barcode{i:02d}": ("dilution_rep1", f"3D7_1e-{i-1}") for i in range(1, 11)},
    **{f"barcode{i:02d}": ("negative_control", "mosquito_only") for i in range(11, 21)},
    **{f"barcode{i:02d}": ("dilution_rep2", f"3D7_1e-{i-21}") for i in range(21, 31)},
    **{f"barcode{i:02d}": ("wild_infected", "field") for i in range(31, 41)},
    **{f"barcode{i:02d}": ("species_id_only", "3D7_spiked_no_plasmo18S_primer") for i in range(41, 51)},
}


def best_alignment_score(bam_path):
    """Return dict: read_name -> (identity, query_coverage) for the best
    (highest identity*coverage) primary+secondary alignment per read."""
    scores = {}
    if not os.path.exists(bam_path):
        return scores
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_supplementary:
                continue
            qlen = read.infer_read_length()
            if not qlen:
                continue
            aligned_len = read.query_alignment_length
            if not aligned_len:
                continue
            nm = read.get_tag("NM") if read.has_tag("NM") else None
            if nm is None:
                continue
            identity = max(0.0, (aligned_len - nm) / aligned_len)
            coverage = aligned_len / qlen
            score = identity * coverage
            name = read.query_name
            if name not in scores or score > scores[name][0] * scores[name][1]:
                scores[name] = (identity, coverage)
    return scores


def classify_sample(barcode, host_suffix="host"):
    host_bam = f"{BAM_DIR}/{barcode}.{host_suffix}.sorted.bam"
    plasmo_bam = f"{BAM_DIR}/{barcode}.plasmo.sorted.bam"
    host_scores = best_alignment_score(host_bam)
    plasmo_scores = best_alignment_score(plasmo_bam)

    all_reads = set(host_scores) | set(plasmo_scores)
    n_plasmo = n_host = n_unclassified = n_total = len(all_reads)

    n_plasmo = 0
    n_host = 0
    n_unclassified = 0
    for name in all_reads:
        h_id, h_cov = host_scores.get(name, (0.0, 0.0))
        p_id, p_cov = plasmo_scores.get(name, (0.0, 0.0))
        h_score = h_id * h_cov
        p_score = p_id * p_cov

        plasmo_ok = p_id >= IDENTITY_THRESH and p_cov >= COVERAGE_THRESH
        host_ok = h_id >= IDENTITY_THRESH and h_cov >= COVERAGE_THRESH

        if plasmo_ok and (not host_ok or p_score > h_score + MARGIN):
            n_plasmo += 1
        elif host_ok and (not plasmo_ok or h_score > p_score + MARGIN):
            n_host += 1
        else:
            n_unclassified += 1

    return {
        "barcode": barcode,
        "n_reads_with_any_alignment": n_total,
        "n_plasmodium": n_plasmo,
        "n_host": n_host,
        "n_unclassified": n_unclassified,
        "pct_plasmodium_of_classified": 100 * n_plasmo / (n_plasmo + n_host) if (n_plasmo + n_host) else float("nan"),
        "pct_plasmodium_of_all_aligned": 100 * n_plasmo / n_total if n_total else float("nan"),
    }


def classify_sample_3way(barcode, host18s_suffix="host18s_standalone"):
    """3-way competitive classification: host score = best of (AgamP4 genome
    alignment, standalone curated/empirical-18S alignment) per read, taken
    separately since merging the curated 18S into the genome index suppressed
    its repetitive k-mers via minimap2's automatic high-occurrence minimizer
    filtering (confirmed: alignment_classification_host18s.tsv was
    numerically identical to the naive genome-only run)."""
    genome_bam = f"{BAM_DIR}/{barcode}.host.sorted.bam"
    standalone18s_bam = f"{BAM_DIR}/{barcode}.{host18s_suffix}.sorted.bam"
    plasmo_bam = f"{BAM_DIR}/{barcode}.plasmo.sorted.bam"

    genome_scores = best_alignment_score(genome_bam)
    standalone18s_scores = best_alignment_score(standalone18s_bam)
    plasmo_scores = best_alignment_score(plasmo_bam)

    all_reads = set(genome_scores) | set(standalone18s_scores) | set(plasmo_scores)

    n_plasmo = n_host = n_unclassified = 0
    n_host_via_18s = 0  # reads rescued specifically by the standalone 18S arm
    for name in all_reads:
        g_id, g_cov = genome_scores.get(name, (0.0, 0.0))
        s_id, s_cov = standalone18s_scores.get(name, (0.0, 0.0))
        p_id, p_cov = plasmo_scores.get(name, (0.0, 0.0))

        g_score = g_id * g_cov
        s_score = s_id * s_cov
        p_score = p_id * p_cov

        # host score = best of the two host arms
        if s_score > g_score:
            h_score, h_id, h_cov, rescued = s_score, s_id, s_cov, True
        else:
            h_score, h_id, h_cov, rescued = g_score, g_id, g_cov, False

        plasmo_ok = p_id >= IDENTITY_THRESH and p_cov >= COVERAGE_THRESH
        host_ok = h_id >= IDENTITY_THRESH and h_cov >= COVERAGE_THRESH

        if plasmo_ok and (not host_ok or p_score > h_score + MARGIN):
            n_plasmo += 1
        elif host_ok and (not plasmo_ok or h_score > p_score + MARGIN):
            n_host += 1
            if rescued and s_score > g_score + MARGIN:
                n_host_via_18s += 1
        else:
            n_unclassified += 1

    n_total = len(all_reads)
    return {
        "barcode": barcode,
        "n_reads_with_any_alignment": n_total,
        "n_plasmodium": n_plasmo,
        "n_host": n_host,
        "n_host_rescued_by_18s_arm": n_host_via_18s,
        "n_unclassified": n_unclassified,
        "pct_plasmodium_of_classified": 100 * n_plasmo / (n_plasmo + n_host) if (n_plasmo + n_host) else float("nan"),
        "pct_plasmodium_of_all_aligned": 100 * n_plasmo / n_total if n_total else float("nan"),
    }


def run_variant_3way(barcodes, out_name, host18s_suffix="host18s_standalone"):
    rows = []
    for bc in barcodes:
        group, label = GROUPS.get(bc, ("unknown", "unknown"))
        res = classify_sample_3way(bc, host18s_suffix=host18s_suffix)
        res["group"] = group
        res["label"] = label
        rows.append(res)
        print(f"{bc:12s} {group:18s} plasmo={res['n_plasmodium']:6d} host={res['n_host']:6d} "
              f"(rescued_by_18s={res['n_host_rescued_by_18s_arm']:6d}) unclass={res['n_unclassified']:6d}  "
              f"%plasmo(of classified)={res['pct_plasmodium_of_classified']:.1f}%")

    fieldnames = ["barcode", "group", "label", "n_reads_with_any_alignment",
                  "n_plasmodium", "n_host", "n_host_rescued_by_18s_arm", "n_unclassified",
                  "pct_plasmodium_of_classified", "pct_plasmodium_of_all_aligned"]
    with open(f"{RESULTS_DIR}/{out_name}", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved {RESULTS_DIR}/{out_name}")
    return rows


def run_variant(barcodes, host_suffix, out_name):
    rows = []
    for bc in barcodes:
        group, label = GROUPS.get(bc, ("unknown", "unknown"))
        res = classify_sample(bc, host_suffix=host_suffix)
        res["group"] = group
        res["label"] = label
        rows.append(res)
        print(f"{bc:12s} {group:18s} plasmo={res['n_plasmodium']:6d} host={res['n_host']:6d} "
              f"unclass={res['n_unclassified']:6d}  %plasmo(of classified)={res['pct_plasmodium_of_classified']:.1f}%")

    fieldnames = ["barcode", "group", "label", "n_reads_with_any_alignment",
                  "n_plasmodium", "n_host", "n_unclassified",
                  "pct_plasmodium_of_classified", "pct_plasmodium_of_all_aligned"]
    with open(f"{RESULTS_DIR}/{out_name}", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved {RESULTS_DIR}/{out_name}")
    return rows


def main():
    with open("barcodes.txt") as f:
        barcodes = [l.strip() for l in f if l.strip()]

    print("=== Variant 1: naive alignment (AgamP4 genome only as host reference) ===")
    run_variant(barcodes, "host", "alignment_classification.tsv")

    print("\n=== Variant 2: improved alignment (AgamP4 genome + curated An. gambiae 18S rRNA, merged index) ===")
    run_variant(barcodes, "host18s", "alignment_classification_host18s.tsv")

    print("\n=== Variant 3: 3-way (AgamP4 genome OR standalone curated-18S index, vs Plasmodium) ===")
    run_variant_3way(barcodes, "alignment_classification_3way.tsv", host18s_suffix="host18s_standalone")

    print("\n=== Variant 4: 3-way (AgamP4 genome OR medaka-polished empirical host-18S, vs Plasmodium) ===")
    run_variant_3way(barcodes, "alignment_classification_4way.tsv", host18s_suffix="host18s_polished")


if __name__ == "__main__":
    main()
