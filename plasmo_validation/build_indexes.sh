#!/bin/bash
# Build the minimap2 indexes used by the Plasmo18S competitive-alignment scripts.
# Commands recovered from logs/index_*.log of the original run.
#   plasmodium_combined.fasta = concatenation of the PlasmoDB-67 P. falciparum 3D7,
#   P. vivax P01, P. malariae UG01 and P. ovale curtisi GH01 genomes.
#   host18s_polished.fasta = Medaka (v2.2.0) polish of host18s_empirical.fasta using
#   mosquito-only negative-control reads (see README).
set -euo pipefail
cd "$(dirname "$0")"

AGAMP4=../VectorBase-67_AgambiaePEST_Genome.fasta
PLASMODB_GENOMES="${PLASMODB_GENOMES:?set PLASMODB_GENOMES to the directory holding the PlasmoDB-67 genome FASTAs}"

cat "$PLASMODB_GENOMES"/PlasmoDB-67_{Pfalciparum3D7,PvivaxP01,PmalariaeUG01,PovalecurtisiGH01}_Genome.fasta > plasmodium_combined.fasta

minimap2 -x map-ont -d agamp4.mmi "$AGAMP4"
minimap2 -x map-ont -d plasmo_combined.mmi plasmodium_combined.fasta
minimap2 -x map-ont -d host18s_standalone.mmi angambiae_18s_all.fasta
minimap2 -x map-ont -d host18s_empirical.mmi host18s_empirical.fasta
minimap2 -x map-ont -d host18s_polished.mmi host18s_polished.fasta
