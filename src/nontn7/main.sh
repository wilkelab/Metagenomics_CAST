#!/usr/bin/env bash
# This script runs the entire non-Tn7 CRISPR-transposon pipeline.

# Halt the pipeline if any errors are encountered
set -euo pipefail

# Directories for input and output
DATA=../../output/nontn7
STORE=$HOME/work
BLASTN_DB="$STORE/trnadb/trnas.fa"
BLASTP_DB="$STORE/uniprotdb/uniprot_combined.fasta"
REBLAST_OUTPUT_DIR="$DATA/reblast"

# Make sure the output directory exists
mkdir -p $DATA
mkdir -p $REBLAST_OUTPUT_DIR

# Find operons that meet our minimal criteria
minimal_file=$DATA/minimal_passing_operons.csv.gz
if [[ ! -e $minimal_file ]]; then
    fd '.*csv' $STORE/missing-effectors $STORE/pipeline-results | parallel -j2 'cat {}' | python fix_paths.py | python rules-all.py | python rules-nocas1-2.py | gzip > $minimal_file
else
    >&2 echo "Minimal file already exists. Skipping..."
fi 

# ===== Class I systems =====
class1_file=$DATA/class1.nuclease_dead.csv.gz
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
#sizes[compact_cas12]="800 3000"
sizes[cas13]="2500 6000"
#sizes[compact_cas13]="800 2500"

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
done | parallel --colsep ' ' "gzip -cd $minimal_file | python rules-{2}.py | python dedup.py | python size_select.py {2} {3} {4} | gzip > $DATA/{1}.csv.gz"

# Make a FASTA file of each effector protein. For Class I proteins, we're using Cas7 at it is essential and does not appear as a fusion with Cas5/6/8.
# For Class 2 proteins, we'll just use the nuclease, either Cas9, Cas12 or Cas13.
# If more than one effector protein is present in an operon, we pick the one with the lowest e-value.

for group in "${!proteins[@]}"; do 
    protein="${proteins[$group]}"
    clustering_file="$DATA/$protein.for_clustering.fasta"
    if [[ ! -e $clustering_file ]]; then
        gzip -cd $DATA/$group.csv.gz | python make_effector_fasta.py $protein > $clustering_file
    fi
done

class1_clustering_file="$DATA/class1.for_clustering.fasta"
if [[ ! -e $class1_clustering_file ]]; then
    gzip -cd $class1_file | python make_effector_fasta.py cas7 > $class1_clustering_file
fi

for group in cas9 cas12 cas13 class1; do
    if [[ ! -e $DATA/$group.clustered_cluster.tsv ]]; then
        tempfile=/tmp/$group.mmseqs.temp
        mmseqs easy-cluster --min-seq-id 0.95 --write-lookup 1 $DATA/$group.for_clustering.fasta $DATA/$group.clustered $tempfile
        rm -r $tempfile
    fi
done















































exit 0  # =====================================================
# Find nuclease-dead Class 2 effectors
for group in "${!proteins[@]}"; do
    fasta="$DATA/$group.effectors.fasta"
    alignment="$DATA/$group.effectors.afa"
    mafft_strategy="$DATA/$group.mafft_strategy"
    protein="${proteins[$group]}"
    residues_file="$DATA/$group.residues.csv"
    nuclease_dead_operons="$DATA/$group.nuclease_dead.csv.gz"
    # Write the protein sequence of every nuclease from each operon to a file
    if [[ ! -e $fasta ]]; then
        exit 1 # THIS IS BROKEN BECAUSE make_effector_fasta isn't the same anymore
        gzip -cd $DATA/$group.csv.gz | python make_effector_fasta.py $protein > $fasta
    fi

    # Perform a multiple sequence alignment so we can identify nucleases with mutated catalytic residues
    if [[ ! -e $alignment ]]; then
        mafft --auto --thread 8 $fasta 2>$alignment.stderr > $alignment
    fi

    # Figure out what the strategy used by MAFFT was
    if [[ ! -e $mafft_strategy ]]; then
        strategy=$(rg --multiline 'Strategy:\n\s(.*)\n' $alignment.stderr --replace '$1') 
        echo "$group $strategy" > $mafft_strategy
        rm $alignment.stderr
    fi

    # Use the alignment and determine which nucleases have (potentially) been inactivated
    if [[ ! -e $residues_file ]]; then
        python identify_catalytic_residues.py ${residues[$protein]} < $alignment > $residues_file 
    fi

    # Gets operons with nuclease-dead effectors and save them to a new file
    if [[ ! -e $nuclease_dead_operons ]]; then
        gzip -cd $DATA/$group.csv.gz | python load_nuclease_dead.py $protein $residues_file | gzip > $nuclease_dead_operons
    fi
done

# Find self-targeting spacers and inverted repeats, both in Class 1 and Class 2 systems
for group in class1 "${!proteins[@]}"; do
    echo $group
done | parallel "gzip -cd $DATA/{}.nuclease_dead.csv.gz | python self_targeting.py | python find_inverted_repeats.py | gzip > $DATA/{}.fully-analyzed.csv.gz"

# Re-BLAST the nuclease-dead systems with Uniref, Swissprot and a tRNA database to see if we've correctly identified the Cas proteins and to understand the genetic context around them
for group in "${!proteins[@]}"; do
    echo $group
done | parallel "gzip -cd $DATA/{}.nuclease_dead.csv.gz | python reblast.py $BLASTN_DB $BLASTP_DB $REBLAST_OUTPUT_DIR/{}"
