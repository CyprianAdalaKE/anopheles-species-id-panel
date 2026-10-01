#!/usr/bin/env bash
# =============================================================================
# 01_compute_coverage.sh (updated for species subfolder structure)
#
# Computes per-locus mean depth for every sample BAM, across all 6 species
# subfolders under aligned/, using your targets.bed file.
# =============================================================================

set -euo pipefail
PANEL_DIR="${PANEL_DIR:-$(cd "$(dirname "$0")" && pwd)}"  # repo/data root; override with env var

# ---- PATHS ----
ALIGNED_DIR="$PANEL_DIR/aligned"
BED_FILE="$PANEL_DIR/nomadsSID.amplicons.bed"
OUT_DIR="$PANEL_DIR/coverage_results"
# ----------------

mkdir -p "$OUT_DIR"

echo "Using targets:"
cat "$BED_FILE"
echo ""

# Loop through each species subfolder
for species_dir in "$ALIGNED_DIR"/*/; do
    species=$(basename "$species_dir")
    echo ">> Species: $species"

    for BAM in "$species_dir"*.sorted.bam; do
        [ -f "$BAM" ] || continue
        sample=$(basename "$BAM" .sorted.bam)

        # Prefix output name with species so samples don't collide across folders
        out_prefix="$OUT_DIR/${species}_${sample}"

        if [ -f "${out_prefix}.regions.bed.gz" ]; then
            echo "   [skip] ${species}_${sample} already processed"
            continue
        fi

        echo "   Processing ${species}_${sample} ..."
        mosdepth \
            --by "$BED_FILE" \
            --no-per-base \
            --thresholds 1,10,30,50,100 \
            "$out_prefix" \
            "$BAM"
    done
done

echo ""
echo "=== Done. ==="
echo "Per-sample coverage files written to: $OUT_DIR"
echo "Next: run 02_aggregate_and_plot.py to build the sample x locus matrix and figures."
