#!/usr/bin/env bash
# Locate each species' MSA-derived amplicon (amplicons_6sp.fa, from extract_amplicons.py) in its own genome;
# the hits_<species>.paf files are the input to insilico_pcr.py.
set -eu
G=${GENOMES:?set GENOMES to the directory with the VectorBase-67 *.k2.fasta genomes}
declare -A GF=([gambiae]=AgambiaePEST [coluzzii]=AcoluzziiAcolN3 [coustani]=AcoustaniAcouGA1
               [arabiensis]=AarabiensisDONGOLA2021 [funestus]=AfunestusAfunGA1 [stephensi]=AstephensiUCISS2018)
for sp in "${!GF[@]}"; do
  awk -v s="|$sp" '/^>/{p=index($0,s)>0} p' amplicons_6sp.fa > tmp_$sp.fa
  minimap2 -x asm20 -t 8 --secondary=no $G/VectorBase-67_${GF[$sp]}_Genome.k2.fasta tmp_$sp.fa 2>/dev/null |
    awk -v OFS='\t' '{print $1,$2,$3,$4,$6,$8,$9,$10,$11,$12}' > hits_$sp.paf
  rm tmp_$sp.fa
done
