#!/usr/bin/env bash
# This script runs the entire non-Tn7 CRISPR-transposon pipeline.

# Halt the pipeline if any errors are encountered
set -euo pipefail

# Directories for input and output
OUTPUT=../../output/nontn7
DATA=../../data
STORE=$HOME/work
BLASTN_DB="$STORE/trnadb/trnas.fa"
BLASTP_DB="$STORE/uniprotdb/uniprot_combined.fasta"
REBLAST_OUTPUT_DIR="$OUTPUT/reblast"

# Set a few parameters
MIN_REBLAST_CLUSTER_SIZE=2
REBLAST_COUNT=3

# Make sure the output directory exists
mkdir -p $OUTPUT
mkdir -p $REBLAST_OUTPUT_DIR

# Find operons that meet our minimal criteria
minimal_file=$OUTPUT/minimal_passing_operons.csv.gz
if [[ ! -e $minimal_file ]]; then
    fd '.*csv' $STORE/missing-effectors $STORE/pipeline-results | parallel -j2 'cat {}' | python fix_paths.py | python rules-all.py | python rules-nocas1-2.py | gzip > $minimal_file
else
    >&2 echo "Minimal file already exists. Skipping..."
fi 

# ===== Class I systems =====
class1_file=$OUTPUT/class1.csv.gz
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
#proteins[compact_cas12]="cas12"
proteins[cas13]="cas13"
#proteins[compact_cas13]="cas13"

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
    if [[ ! -e "$OUTPUT/$group.csv.gz" ]]; then 
        protein_name="${proteins[$group]}"
        protein_range="${sizes[$group]}"
        echo $group $protein_name $protein_range
    else
        >&2 echo "$group operons already exist. Skipping..."
    fi
done | parallel --colsep ' ' "gzip -cd $minimal_file | python rules-{2}.py | python dedup.py | python size_select.py {2} {3} {4} | gzip > $OUTPUT/{1}.csv.gz"

# Make a FASTA file of each effector protein. For Class I proteins, we're using Cas7 at it is essential and does not appear as a fusion with Cas5/6/8.
# For Class 2 proteins, we'll just use the nuclease, either Cas9, Cas12 or Cas13.
# If more than one effector protein is present in an operon, we pick the one with the lowest e-value.
# We also make a separate file for the Class 2 systems, which has a reference protein at the top of the file. We will perform a multiple sequence alignment 
# on this file to determine which proteins are nuclease-active.

for group in "${!proteins[@]}"; do 
    protein="${proteins[$group]}"
    clustering_file="$OUTPUT/$group.for_clustering.fasta"
    if [[ ! -e $clustering_file ]]; then
        # Make the FASTA file to be used for clustering
        gzip -cd $OUTPUT/$group.csv.gz | python make_effector_fasta.py $protein > $clustering_file
        # Make the FASTA file to be used for alignment
        cat $DATA/nontn7/$group.fasta $clustering_file > $OUTPUT/$group.for_alignment.fasta
    fi
done

class1_clustering_file="$OUTPUT/class1.for_clustering.fasta"
if [[ ! -e $class1_clustering_file ]]; then
    gzip -cd $class1_file | python make_effector_fasta.py cas7 > $class1_clustering_file
fi


for group in "${!proteins[@]}"; do
    fasta=$OUTPUT/$group.for_alignment.fasta
    alignment=$OUTPUT/$group.afa
    mafft_strategy="$OUTPUT/$group.mafft_strategy"
    protein="${proteins[$group]}"
    residues_file="$OUTPUT/$group.residues.csv"
    nuclease_dead_operons="$OUTPUT/$group.nuclease_dead.csv.gz"
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
        gzip -cd $OUTPUT/$group.csv.gz | python load_nuclease_dead.py $protein $residues_file | gzip > $nuclease_dead_operons
    fi
done


for group in cas9 cas12 cas13 class1; do
    all_seqs_file="$OUTPUT/$group.clustered_all_seqs.fasta"
    clustering_file="$OUTPUT/$group.for_clustering.fasta"
    if [[ ! -e $all_seqs_file ]]; then
        tempfile=/tmp/$group.mmseqs.temp
        mmseqs easy-cluster --min-seq-id 0.95 --write-lookup 1 $clustering_file $OUTPUT/$group.clustered $tempfile
        rm -r $tempfile
        mkdir -p $OUTPUT/$group.clusters
        gzip -cd $OUTPUT/$group.csv.gz | python get_clustered_operons.py $all_seqs_file $OUTPUT/$group.clusters
    fi
done


for group in "${!proteins[@]}"; do
    ndcluster_dir=$OUTPUT/$group.nuclease_dead_clusters
    if [[ -e $ndcluster_dir ]]; then
        echo "nuclease-dead analysis done for $group"
        continue
    fi
    mkdir -p $ndcluster_dir
    # using the nuclease-active list from above, determine which clusters have no nuclease-active proteins. write that to a machine-readable list
    nuclease_dead_operons="$OUTPUT/$group.nuclease_dead.csv.gz"
    cluster_gzip_dir=$OUTPUT/$group.clusters
    for nuclease_dead_cluster_filename in $(python count_nuclease_dead_in_cluster.py $nuclease_dead_operons $cluster_gzip_dir); do
        cp $cluster_gzip_dir/$nuclease_dead_cluster_filename $ndcluster_dir
    done
done


for group in "${!proteins[@]}"; do
    for filename in $(ls $OUTPUT/$group.nuclease_dead_clusters | rg '.*csv.gz'); do
        plot_output_dir=$OUTPUT/$group.nuclease_dead_clusters/plots/$filename
        if [[ -e $plot_output_dir ]]; then
            continue
        fi
        mkdir -p $plot_output_dir
        gzip -cd $OUTPUT/$group.nuclease_dead_clusters/$filename | python plot_operons.py $plot_output_dir
    done
done

 # Find self-targeting spacers and inverted repeats in Class 2 systems
for group in "${!proteins[@]}"; do
    fully_analyzed_dir=$OUTPUT/$group.fully-analyzed
    if [[ -e $fully_analyzed_dir ]]; then
        >&2 echo "Inverted repeats and self-targeting spacers already found for $group"
        continue
    fi
    mkdir -p $fully_analyzed_dir
    for filename in $(ls $OUTPUT/$group.nuclease_dead_clusters | rg '.*csv.gz'); do
        echo $group $filename
    done
done | parallel -j 8 --colsep ' ' "gzip -cd $OUTPUT/{1}.nuclease_dead_clusters/{2} | python self_targeting.py | python find_inverted_repeats.py | gzip > $OUTPUT/{1}.fully-analyzed/{2}"

# Find self-targeting spacers and inverted repeats in Class 1 systems
fully_analyzed_dir=$OUTPUT/class1.fully-analyzed
if [[ ! -e $fully_analyzed_dir ]]; then
    mkdir -p $fully_analyzed_dir
    for filename in $(ls $OUTPUT/class1.clusters | rg '.*csv.gz'); do
        echo $filename
    done | parallel -j 8 --colsep ' ' "gzip -cd $OUTPUT/class1.clusters/{1} | python self_targeting.py | python find_inverted_repeats.py | gzip > $fully_analyzed_dir/{1}"
else
    >&2 echo "Inverted repeats and self-targeting spacers already found for Class 1"
fi


exit 99

#for group in "${!proteins[@]}"; do
for group in cas9 class1; do
    filedir=$OUTPUT/$group.fully-analyzed
    for filename in $(ls $filedir); do
        cluster_directory=$(basename $filename)
        directory="$REBLAST_OUTPUT_DIR/$group/$cluster_directory"
        mkdir -p $directory
        completed_file=$REBLAST_OUTPUT_DIR/$group/completed.txt
        touch $completed_file

        if [[ ! $(rg $filename $completed_file) ]]; then
            >&2 echo "Reblasting $group operons from $filename"
            # Perform the re-BLASTing. If we skip the file because it has too few clusters, we print nothing so that we can lower the minimum cluster size later. If we do print the filename it indicates that we BLASTed as many operons as we wanted.
            # TODO: rename newreblast.py to reblast.py once Class 1 is done BLASTing
            python newreblast.py $BLASTN_DB $BLASTP_DB $MIN_REBLAST_CLUSTER_SIZE $REBLAST_COUNT $directory $OUTPUT/$group.clusters/$filename >> $completed_file
        else
            >&2 echo "Already finished $group operons from $filename"
        fi
    done
done
