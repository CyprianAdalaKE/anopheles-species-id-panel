#!/usr/bin/env bash
# =============================================================================
# run_kraken_bracken.sh
#
# Runs Kraken2 taxonomic classification followed by Bracken abundance
# re-estimation on every sample across all 6 Anopheles species folders and
# the Plasmodium infectivity folder, using the pre-built custom database.
# =============================================================================

set -euo pipefail
PANEL_DIR="${PANEL_DIR:-$(cd "$(dirname "$0")" && pwd)}"  # repo/data root; override with env var

# ---- PATHS ----
BASE="$PANEL_DIR"
DB="${KRAKEN_DB:?set KRAKEN_DB to the Kraken2 database directory}"
OUT_DIR="$BASE/kraken_results"
THREADS=4
READ_LEN=100
BRACKEN_LEVEL="S"
BRACKEN_BIN="${BRACKEN_BIN:-bracken}"
# ----------------

# Map each source folder to a clean label (same as alignment stage)
declare -A SOURCE_DIRS=(
    ["gam_final"]="gambiae"
    ["arabiensis"]="arabiensis"
    ["colu_final"]="coluzzii"
    ["fun_final"]="funestus"
    ["coustani_combined"]="coustani"
    ["steph_combined"]="stephensi"
    ["infectivity"]="infectivity"
)

mkdir -p "$OUT_DIR"

for dir in "${!SOURCE_DIRS[@]}"; do
    label="${SOURCE_DIRS[$dir]}"
    input_dir="$BASE/$dir"
    out_dir="$OUT_DIR/$label"
    mkdir -p "$out_dir"

    echo ""
    echo ">>> Processing: $label (folder: $dir)"

    shopt -s nullglob
    files=("$input_dir"/*barcode*)
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        echo "    WARNING: no barcode files found in $input_dir -- skipping."
        continue
    fi

    for f in "${files[@]}"; do
        fname=$(basename "$f")
        sample=$(echo "$fname" | sed -e 's/^combined\.//' -e 's/\.fastq[0-9]*$//')

        kreport="$out_dir/${sample}.kreport"
        bracken_out="$out_dir/${sample}.bracken"

        if [ -f "$bracken_out" ]; then
            echo "    [skip] $label/$sample already classified"
            continue
        fi

        echo "    Classifying $label/$sample ..."

        if [ -f "$kreport" ]; then
            echo "        (kraken2 report already exists, skipping kraken2)"
        else
            kraken2 \
                --db "$DB" \
                --threads "$THREADS" \
                --report "$kreport" \
                --output "$out_dir/${sample}.kraken2" \
                "$f" > /dev/null
        fi

        # Bracken re-estimates abundance at species level using the kraken report
        "$BRACKEN_BIN" \
            -d "$DB" \
            -i "$kreport" \
            -o "$bracken_out" \
            -w "$out_dir/${sample}_bracken.kreport" \
            -r "$READ_LEN" \
            -l "$BRACKEN_LEVEL" \
            -t 10
    done
done

echo ""
echo "=== Done. ==="
echo "Kraken2/Bracken results organized under: $OUT_DIR/<label>/"
echo "Each sample has: <sample>.kreport, <sample>.kraken2, <sample>.bracken, <sample>_bracken.kreport"
