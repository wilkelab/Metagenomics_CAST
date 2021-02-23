#!/usr/bin/env bash

# Aligns all Cas12k proteins from NCBI (plus two nuclease-active control Cas12a proteins) and confirms that our nuclease-dead analysis works.

set -euo pipefail

OUTPUT=../../output/nontn7
DATA=../../data/nontn7

alignment=$OUTPUT/cas12k.afa
residues_file=$OUTPUT/cas12k.residues.csv
mafft --auto --thread 8 $DATA/ncbi-cas12k-seqs.fasta 2>/dev/null > $alignment
python identify_catalytic_residues.py D908 E993 < $alignment > $residues_file 
