#!/usr/bin/env bash
# This script runs the Tn7-class2 CRISPR-transposon pipeline.

# Halt the pipeline if any errors are encountered
set -euo pipefail

# Handle command line arguments
data_partition=$1
if [[ "$data_partition" == "ncbi" ]]; then
    DATA="/stor/work/Wilke/crisprtn7/NCBI_gene_finder.csv"
elif [[ "$data_partition" == "meta" ]]; then
    DATA="/stor/work/Wilke/crisprtn7/typeI_Meta_gene_finder.csv"
else
    echo "Data partition doesn't exist. Exiting..." && exit
fi

# Output directory
#current_time=$(date "+%Y%m%d-%H")
RESULTS="./results"

# Make sure the output directory exists
mkdir -p $RESULTS


Operon_analyzer_file_cas9=$RESULTS/tn7_cas9_$data_partition.csv
if [[ ! -e $Operon_analyzer_file_cas9 ]]; then
    python rules-cas9.py $DATA $Operon_analyzer_file_cas9
else
    >&2 echo "Tn7_cas9 file already exists. Skipping..."
fi

# Find candidates that meet the cas12 criteria
Operon_analyzer_file_cas12=$RESULTS/tn7_cas12_$data_partition.csv
if [[ ! -e $Operon_analyzer_file_cas12 ]]; then
    python rules-cas12.py $DATA $Operon_analyzer_file_cas12
else
    >&2 echo "Tn7_cas12 file already exists. Skipping..."
fi

# Find candidates that meet the cas13 criteria
Operon_analyzer_file_cas13=$RESULTS/tn7_cas13_$data_partition.csv
if [[ ! -e $Operon_analyzer_file_cas13 ]]; then
    python rules-cas13.py $DATA $Operon_analyzer_file_cas13
else
    >&2 echo "Tn7_cas13 file already exists. Skipping..."
fi

# Find candidates that meet the cas14 criteria
Operon_analyzer_file_cas14=$RESULTS/tn7_cas14_$data_partition.csv
if [[ ! -e $Operon_analyzer_file_cas14 ]]; then
    python rules-cas14.py $DATA $Operon_analyzer_file_cas14
else
    >&2 echo "Tn7_cas14 file already exists. Skipping..."
fi

# ===== Operon Visualization =====
# For testing/debugging - comment out as needed
diagram_dir=$RESULTS/class2_diagrams_$data_partition
if [[ -e $diagram_dir ]]; then
    rm -rf $diagram_dir/*
    mkdir $diagram_dir/tn7_cas9
    mkdir $diagram_dir/tn7_cas12
    mkdir $diagram_dir/tn7_cas13
    mkdir $diagram_dir/tn7_cas14
else
    mkdir $diagram_dir
    mkdir $diagram_dir/tn7_cas9
    mkdir $diagram_dir/tn7_cas12
    mkdir $diagram_dir/tn7_cas13
    mkdir $diagram_dir/tn7_cas14
fi

python png_generator.py $Operon_analyzer_file_cas9 $DATA $diagram_dir/tn7_cas9
python png_generator.py $Operon_analyzer_file_cas12 $DATA $diagram_dir/tn7_cas12
python png_generator.py $Operon_analyzer_file_cas13 $DATA $diagram_dir/tn7_cas13
python png_generator.py $Operon_analyzer_file_cas14 $DATA $diagram_dir/tn7_cas14
