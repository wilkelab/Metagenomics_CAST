#!/usr/bin/env bash

# This script calculates several numerical values that we report in the text of the paper.

# TODO: this is no longer valid since we've split things up into ±CRISPR arrays, and also all/composite/tn3
OUTPUT="/stor/home/rybarskj/Metagenomics_CRISPR_transposons/output/nontn7"
GROUPS="cas9 cas12 cas13 class1"
# 1. The count of the number of systems we found that meet the basic rules.

function count_operons () {
    echo $(fd -a '.*csv.gz' $OUTPUT/$1.clusters | parallel -j2 'gzip -cd {}' | python count_operons.py)
}

cas9count=$(count_operons cas9)
cas12count=$(count_operons cas12)
cas13count=$(count_operons cas13)
class1count=$(count_operons class1)

printf "The pipeline identified %'d putative Cas9-encoding systems, %'d Cas12-encoding systems, %'d Cas13-encoding systems and %'d Type I systems.\n" $cas9count $cas12count $cas13count $class1count 
