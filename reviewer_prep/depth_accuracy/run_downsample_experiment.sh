#!/bin/bash
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate anopheles-panel
cd "$(dirname "$0")"

DB="${KRAKEN_DB:?set KRAKEN_DB to the Kraken2 database directory}"
DEPTHS="50 100 200 500 1000 2000 5000"
BRACKEN="${BRACKEN_BIN:-bracken}"

mkdir -p subsampled_fastq kraken_out
> depth_experiment_results.tsv
echo -e "species\tbarcode\ttarget_depth\tactual_reads\tpredicted_species\tcorrect" > depth_experiment_results.tsv

run_one() {
    sp="$1"; bc="$2"; native="$3"; fq="../../$4"; depth="$5"
    tag="${sp}_${bc}_d${depth}"
    sub_fq="subsampled_fastq/${tag}.fastq"
    kreport="kraken_out/${tag}.kreport"
    kraken2out="kraken_out/${tag}.kraken2"
    bracken_out="kraken_out/${tag}.bracken"

    if [ "$depth" -ge "$native" ]; then
        actual="$native"
        cp "$fq" "$sub_fq"
    else
        actual="$depth"
        seqtk sample -s100 "$fq" "$depth" > "$sub_fq" 2>>kraken_out/seqtk.log
        if [ ! -s "$sub_fq" ]; then
            echo "EMPTY OUTPUT: seqtk sample -s100 $fq $depth (tag=$tag)" >> kraken_out/seqtk.log
        fi
    fi

    kraken2 --db "$DB" --threads 2 --report "$kreport" --output "$kraken2out" "$sub_fq" > /dev/null 2>>kraken_out/kraken.log
    "$BRACKEN" -d "$DB" -i "$kreport" -o "$bracken_out" -r 100 -l S > /dev/null 2>>kraken_out/bracken.log

    pred="NA"
    if [ -s "$bracken_out" ]; then
        pred=$(tail -n +2 "$bracken_out" | sort -t$'\t' -k6,6 -rn | head -1 | cut -f1)
    fi
    pred_lc=$(echo "$pred" | tr '[:upper:]' '[:lower:]')
    correct="FALSE"
    case "$pred_lc" in
        *"$sp"*) correct="TRUE" ;;
    esac

    echo -e "${sp}\t${bc}\t${depth}\t${actual}\t${pred}\t${correct}" >> depth_experiment_results.tsv
    rm -f "$sub_fq"
}
export -f run_one
export DB BRACKEN

tail -n +2 selected_samples.tsv | while IFS=$'\t' read -r sp bc native fq; do
    for d in $DEPTHS; do
        echo "$sp $bc $native $fq $d"
    done
done | xargs -P 4 -L 1 bash -c 'run_one "$@"' _

echo "ALL DONE" >> depth_experiment_results.tsv
