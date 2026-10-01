#!/usr/bin/env bash
# =============================================================================
# build_database.sh
#
# Builds the custom Kraken2 + Bracken database (species_id_infectivity_db)
# used for species identification and Plasmodium detection.
#
# Reference genomes (VectorBase / PlasmoDB release 67) must be downloaded into
# $DB/genomes/ first. Because they are not from NCBI RefSeq, each FASTA header
# is tagged with kraken:taxid|<taxid> by add_kraken_to_fasta.py.
#
# NCBI taxonomy (names.dmp, nodes.dmp, nucl_gb.accession2taxid) must be present
# in $DB/taxonomy/ (kraken2-build --download-taxonomy --db $DB).
# =============================================================================

set -euo pipefail

# ---- PATHS ----
DB="${KRAKEN_DB:?set KRAKEN_DB to the database directory}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
THREADS=8
READ_LEN=100
# ----------------

# genome FASTA basename -> NCBI taxid
declare -A TAXIDS=(
    ["VectorBase-67_AgambiaePEST_Genome"]=7165
    ["VectorBase-67_AcoluzziiAcolN3_Genome"]=1518534
    ["VectorBase-67_AarabiensisDONGOLA2021_Genome"]=7173
    ["VectorBase-67_AfunestusAfunGA1_Genome"]=62324
    ["VectorBase-67_AstephensiUCISS2018_Genome"]=30069
    ["VectorBase-67_AcoustaniAcouGA1_Genome"]=139045
    ["PlasmoDB-67_Pfalciparum3D7_Genome"]=5833
    ["PlasmoDB-67_PvivaxP01_Genome"]=5855
    ["PlasmoDB-67_PmalariaeUG01_Genome"]=5863
    ["PlasmoDB-67_PovalecurtisiGH01_Genome"]=5800
)

for name in "${!TAXIDS[@]}"; do
    python "$SCRIPT_DIR/add_kraken_to_fasta.py" "$DB/genomes/${name}.fasta" "${TAXIDS[$name]}" "$DB/genomes" > /dev/null
    kraken2-build --add-to-library "$DB/genomes/${name}.k2.fasta" --db "$DB"
done

# Default parameters: k = 35, minimizer length = 31
kraken2-build --build --db "$DB" --threads "$THREADS"

# Bracken k-mer distribution for 100 bp reads
bracken-build -d "$DB" -t "$THREADS" -k 35 -l "$READ_LEN"
