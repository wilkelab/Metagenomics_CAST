#!/bin/bash

# Counts how many operons are in each gzipped file in a directory.

for filename in $(fd csv.gz ../../output/nontn7/$1); do 
    # Find all unique combinations of accession and contig coordinates
    count=$(gzip -cd $filename | awk 'BEGIN {FS=","};{print $1 $2}' | uniq | wc -l)
    echo "$count,$filename"
done
