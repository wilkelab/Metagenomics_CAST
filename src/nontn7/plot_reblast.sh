#!/usr/bin/env bash


set -euo pipefail

OUTPUT=/stor/home/rybarskj/Metagenomics_CRISPR_transposons/output/nontn7

for group in class1; do
    echo "Group: $group"
    fd '.*csv$' $OUTPUT/reblast/$group | parallel -j2 'cat {}' | python plot_reblast.py $OUTPUT/$group.fully-analyzed $OUTPUT/plots/$group
done
