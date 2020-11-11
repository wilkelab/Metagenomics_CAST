#!/usr/bin/env bash
# This script runs the Tn7-class1 CRISPR-transposon pipeline.

# Halt the pipeline if any errors are encountered
set -euo pipefail

# Handle command line arguments
data_partition=$1
if [[ "$data_partition" == "ncbi" ]]; then
    DATA="/stor/scratch/Wilke/amh7958/crispr-transposons/data/ncbi"
elif [[ "$data_partition" == "meta" ]]; then
    DATA="/stor/scratch/Wilke/amh7958/crispr-transposons/data/metagenomic"
else
    echo "Data partition doesn't exist. Exiting..." && exit
fi

# Output directory
#current_time=$(date "+%Y%m%d-%H")
RESULTS="/stor/scratch/Wilke/amh7958/crispr-transposons/results"

# Make sure the output directory exists
mkdir -p $RESULTS

# Find candidates that meet the minimal criteria
minimal_file=$RESULTS/minimal_tn7_$data_partition.csv
if [[ ! -e $minimal_file ]]; then
    find $DATA -name "*.csv" | parallel -j2 'cat {}' | python3 rules-minimal.py > $minimal_file
else
    >&2 echo "Minimal file already exists. Skipping..."
fi

# Apply rules to select only good-quality Tn7-class1 loci
class1_file=$RESULTS/class1_tn7_$data_partition.csv
if [[ ! -e $class1_file ]]; then
    cat $minimal_file | python3 rules-tn7-class1.py | python3 ../nontn7/dedup.py > $class1_file 
else
    >&2 echo "Class I operons already exist. Skipping..."
fi

# ===== Operon Visualization =====
# For testing/debugging - comment out as needed
diagram_dir=$RESULTS/class1_diagrams_$data_partition
if [[ -e $diagram_dir ]]; then
    rm $diagram_dir/*.png
else
    mkdir $diagram_dir
fi

cat $class1_file | python3 operon-diagrams.py $diagram_dir
