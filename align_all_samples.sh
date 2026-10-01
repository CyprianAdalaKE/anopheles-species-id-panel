#!/usr/bin/env bash
# =============================================================================
# align_all_samples.sh
#
# Aligns every FASTQ sample across all 6 Anopheles species folders to the
# AgamP4 reference genome using minimap2, producing sorted, indexed BAM
# files organized by species.
#
# Handles 3 different filename patterns automatically:
#   - barcode01.fastq ... barcode80.fastq         (gam_final, arabiensis, colu_final, fun_final)
#   - combined.barcode01.fastq ... barcode80.fastq (coustani_combined)
#   - combined.barcode01 ... barcode80 (no ext)    (steph_combined)
# =============================================================================

set -euo pipefail
PANEL_DIR="${PANEL_DIR:-$(cd "$(dirname "$0")" && pwd)}"  # repo/data root; override with env var

# ---- PATHS (already confirmed) ----
BASE="$PANEL_DIR"
REF="$BASE/VectorBase-67_AgambiaePEST_Genome.fasta"
OUT_BASE="$BASE/aligned"
THREADS=4   # increase if your machine has more cores (check with: nproc)
# ------------------------------------

# Map each input folder to a clean species label for output organization
declare -A SPECIES_DIRS=(
    ["gam_final"]="gambiae"
    ["arabiensis"]="arabiensis"
    ["colu_final"]="coluzzii"
    ["fun_final"]="funestus"
    ["coustani_combined"]="coustani"
    ["steph_combined"]="stephensi"
)

echo "=== Step 1: Indexing reference genome (only needs to happen once) ==="
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
    echo "Reference indexed."
else
    echo "Reference already indexed, skipping."
fi

mkdir -p "$OUT_BASE"

echo ""
echo "=== Step 2: Aligning samples ==="

for dir in "${!SPECIES_DIRS[@]}"; do
    species="${SPECIES_DIRS[$dir]}"
    input_dir="$BASE/$dir"
    out_dir="$OUT_BASE/$species"
    mkdir -p "$out_dir"

    echo ""
    echo ">>> Species: $species (folder: $dir)"

    # Only pick up files that look like barcode samples, regardless of
    # naming pattern or extension -- avoids accidentally grabbing stray
    # .tsv or other files sitting in the same folder.
    shopt -s nullglob
    files=("$input_dir"/*barcode*)
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        echo "    WARNING: no barcode files found in $input_dir -- skipping."
        continue
    fi

    for f in "${files[@]}"; do
        fname=$(basename "$f")

        # Normalize sample name: strip leading "combined." and trailing ".fastq" if present
        sample=$(echo "$fname" | sed -e 's/^combined\.//' -e 's/\.fastq$//')

        out_bam="$out_dir/${sample}.sorted.bam"

        if [ -f "$out_bam" ]; then
            echo "    [skip] $species/$sample already aligned"
            continue
        fi

        echo "    Aligning $species/$sample ..."
        minimap2 -ax map-ont -t "$THREADS" "$REF" "$f" 2>/dev/null \
            | samtools sort -@ "$THREADS" -o "$out_bam" -
        samtools index "$out_bam"
    done
done

echo ""
echo "=== Done. ==="
echo "Aligned BAMs organized under: $OUT_BASE/<species>/<sample>.sorted.bam"
echo "Next step: per-locus coverage analysis on these BAMs."
