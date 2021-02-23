#!/usr/bin/env bash

# This script calculates several numerical values that we report in the text of the paper.

OUTPUT="/stor/home/rybarskj/Metagenomics_CRISPR_transposons/output/nontn7/systems-with-crispr-arrays/all"
GROUPS="cas9 cas12 cas13 class1"

function count_operons () {
    echo $(fd -a '.*csv.gz' $OUTPUT/$1.clusters | parallel -j2 'gzip -cd {}' | python count_operons.py)
}

function count_nd_operons() {
    echo $(fd -a '.*csv.gz' $OUTPUT/$1.fully-analyzed | parallel -j2 'gzip -cd {}' | python count_operons.py)
}

# 1. The count of the number of systems we found that meet the basic rules.

cas9count=$(count_operons cas9)
cas12count=$(count_operons cas12)
cas13count=$(count_operons cas13)
class1count=$(count_operons class1)

printf "The pipeline identified %'d putative Cas9-encoding systems, %'d Cas12-encoding systems, %'d Cas13-encoding systems and %'d Type I systems.\n" $cas9count $cas12count $cas13count $class1count 

deduplicated_count=$(gzip -cd /stor/home/rybarskj/Metagenomics_CRISPR_transposons/output/nontn7/minimal_passing_systems_with_crispr_arrays.csv.gz | python dedup.py | python count_operons.py)
original_count=$(gzip -cd /stor/home/rybarskj/Metagenomics_CRISPR_transposons/output/nontn7/minimal_passing_systems_with_crispr_arrays.csv.gz | python count_operons.py)
printf "There were %'d systems that passed our minimal rules, and %'d systems after deduplication.\n" $original_count $deduplicated_count


# 2. Determine the percentage of systems that were nuclease-dead.

cas9nd_count=$(count_nd_operons cas9)
cas9_nd_ratio=$(echo "scale=4; $cas9nd_count/$cas9count*100" | bc)
cas12nd_count=$(count_nd_operons cas12)
cas12_nd_ratio=$(echo "scale=4; $cas12nd_count/$cas12count*100" | bc)
cas13nd_count=$(count_nd_operons cas13)
cas13_nd_ratio=$(echo "scale=4; $cas13nd_count/$cas13count*100" | bc)

printf "We found %f%% of Cas9, %f%% of Cas12, and %f%% of Cas13 had nuclease-inactivating mutations or deletions.\n" $cas9_nd_ratio $cas12_nd_ratio $cas13_nd_ratio
