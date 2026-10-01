#!/usr/bin/env bash
# =============================================================================
# align_infectivity.sh
#
# Aligns the 80 "infectivity" samples (Plasmodium 18S rRNA amplicon target)
# to the P. falciparum 3D7 reference genome -- kept completely separate from
# the Anopheles species BAMs, since this is a different target and reference.
# =============================================================================

set -euo pipefail
PANEL_DIR="${PANEL_DIR:-$(cd "$(dirname "$0")" && pwd)}"  # repo/data root; override with env var

# ---- PATHS ----
INPUT_DIR="$PANEL_DIR/infectivity"
REF="${PF_REF:?set PF_REF to PlasmoDB-67_Pfalciparum3D7_Genome.fasta}"
OUT_DIR="$PANEL_DIR/aligned_plasmodium"
THREADS=4
# ----------------

echo "=== Step 1: Checking reference index ==="
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
    echo "Reference indexed."
else
    echo "Reference already indexed, skipping."
fi

mkdir -p "$OUT_DIR"

echo ""
echo "=== Step 2: Aligning infectivity samples to Pf3D7 ==="

shopt -s nullglob
files=("$INPUT_DIR"/*barcode*)
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: no barcode files found in $INPUT_DIR"
    exit 1
fi

for f in "${files[@]}"; do
    fname=$(basename "$f")

    # Normalize sample name: strip leading "combined." and trailing ".fastq" if present
    sample=$(echo "$fname" | sed -e 's/^combined\.//' -e 's/\.fastq$//')

    out_bam="$OUT_DIR/${sample}.sorted.bam"

    if [ -f "$out_bam" ]; then
        echo "[skip] $sample already aligned"
        continue
    fi

    echo "Aligning $sample ..."
    minimap2 -ax map-ont -t "$THREADS" "$REF" "$f" 2>/dev/null \
        | samtools sort -@ "$THREADS" -o "$out_bam" -
    samtools index "$out_bam"
done

echo ""
echo "=== Done. ==="
echo "Plasmodium-aligned BAMs written to: $OUT_DIR"
