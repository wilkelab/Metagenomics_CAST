#!/bin/bash

# Counts how many operons passed all the rules for each branch of rules

for filename in $(fd csv.gz ../../output/nontn7); do 
    # Find all unique combinations of accession and contig coordinates
    count=$(gzip -cd $filename | awk 'BEGIN {FS=","};{print $1 $2}' | uniq | wc -l)
    echo "$count,$filename"
done
