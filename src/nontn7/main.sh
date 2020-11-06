#!/usr/bin/env bash
# This script runs the entire non-Tn7 CRISPR-transposon pipeline.

# Halt the pipeline if any errors are encountered
set -euo pipefail

# Directories for input and output
DATA=data
RESULTS=$HOME/work
BLASTN_DB="/stor/work/Wilke/rybarskj/trnadb/trnas.fa"
BLASTP_DB="/stor/work/Wilke/rybarskj/uniprotdb/uniprot_combined.fasta"
REBLAST_OUTPUT_DIR="$DATA/reblast"

# Make sure the output directory exists
mkdir -p $DATA
mkdir -p $REBLAST_OUTPUT_DIR

# Find operons that meet our minimal criteria
minimal_file=$DATA/minimal_passing_operons.csv.gz
if [[ ! -e $minimal_file ]]; then
    fd '.*csv' $RESULTS/missing-effectors $RESULTS/pipeline-results | parallel -j2 'cat {}' | python fix-paths.py | python rules-all.py | python rules-nocas1-2.py | gzip > $minimal_file
else
    >&2 echo "Minimal file already exists. Skipping..."
fi 


# ===== Class I systems =====
class1_file=$DATA/class1.csv.gz
if [[ ! -e $class1_file ]]; then
    gzip -cd $minimal_file | python rules-nocas3-10.py | python rules-class1.py | python dedup.py | gzip > $class1_file 
else
    >&2 echo "Class I operons already exist. Skipping..."
fi


# ===== Class II systems =====
# We split certain systems by size but some scripts need the exact name of the protein
declare -A proteins
proteins[cas9]="cas9"
proteins[cas12]="cas12"
proteins[compact_cas12]="cas12"
proteins[cas13]="cas13"
proteins[compact_cas13]="cas13"

# Sizes are the minimum and maximum length of the effector gene, in base pairs
declare -A sizes
sizes[cas9]="2000 6000"
sizes[cas12]="3000 6000"
sizes[compact_cas12]="800 3000"
sizes[cas13]="2500 6000"
sizes[compact_cas13]="800 2500"

# These are the catalytic residues in the reference proteins
declare -A residues
residues[cas9]="D10 H840"
residues[cas12]="D908 E993"
residues[cas13]="R472 H477 R1048 H1053"

# Apply type-specific rules, deduplicate operons, and select systems with appropriately-sized effectors
for group in "${!sizes[@]}"; do
    if [[ ! -e "$DATA/$group.csv.gz" ]]; then 
        protein_name="${proteins[$group]}"
        protein_range="${sizes[$group]}"
        echo "$group $protein_name $protein_range"
    else
        >&2 echo "$group operons already exist. Skipping..."
    fi
done | parallel --colsep ' ' "gzip -cd $minimal_file | python rules-{2}.py | python dedup.py | python size-select.py {2} {3} {4} | gzip > $DATA/{1}.csv.gz"

# Select operons that may have nuclease-dead effectors
for group in "${!proteins[@]}"; do
    fasta="$DATA/$group.effectors.fasta"
    alignment="$DATA/$group.effectors.afa"
    protein="${proteins[$group]}"
    residues_file="$DATA/$group.residues.csv"
    nuclease_dead_operons="$DATA/$group.nuclease_dead.csv.gz"
    # Write the protein sequence of every nuclease from each operon to a file
    if [[ ! -e $fasta ]]; then
        gzip -cd $DATA/$group.csv.gz | python make-effector-fasta.py $protein > $fasta
    fi

    # Perform a multiple sequence alignment so we can identify nucleases with mutated catalytic residues
    if [[ ! -e $alignment ]]; then
        mafft --auto --thread 8 $fasta 2>/dev/null > $alignment
    fi

    # Use the alignment and determine which nucleases have (potentially) been inactivated
    if [[ ! -e $residues_file ]]; then
        python identify-catalytic-residues.py ${residues[$protein]} < $alignment > $residues_file 
    fi

    if [[ ! -e $nuclease_dead_operons ]]; then
        gzip -cd $DATA/$group.csv.gz | python load-nuclease-dead.py $protein $residues_file | gzip > $nuclease_dead_operons
    fi
done

for group in "${!proteins[@]}"; do
    echo $group
done | parallel "gzip -cd $DATA/{}.nuclease_dead.csv.gz | python reblast.py $BLASTN_DB $BLASTP_DB $REBLAST_OUTPUT_DIR/{}"
