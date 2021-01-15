#!/usr/bin/env bash
# This script runs the entire non-Tn7 CRISPR-transposon pipeline.


set -euo pipefail  # Halt the pipeline if any errors are encountered

OUTPUT=../../output/nontn7
DATA=../../data
STORE=$HOME/work
INPUT=$(ls /stor/work/Wilke/amh7958/pipeline-results/*tar.gz /stor/work/Wilke/amh7958/pipeline-results/missing_effectors/*.tar.gz)
BLASTN_DB="$STORE/trnadb/trnas.fa"
BLASTP_DB="$STORE/uniprotdb/uniprot_combined.fasta"
REBLAST_OUTPUT_DIR="$OUTPUT/reblast"
MIN_REBLAST_CLUSTER_SIZE=1
REBLAST_COUNT=3
KEEP_PATHS="NO"  # change to "YES" if you're running this pipeline on some other system than the Wilke cluster. You will need to update the values for $STORE and $INPUT, and probably most of the hard-coded values listed above.

# TODO: This looks ridiculous because we removed the weird corner cases it allowed us to use.
# It can be replaced by a simple list of strings probably.
declare -A proteins
proteins[cas9]="cas9"
proteins[cas12]="cas12"
proteins[cas13]="cas13"

# Sizes are the minimum and maximum length of the effector genes that we allow, in base pairs
declare -A sizes
sizes[cas9]="2000 6000"
sizes[cas12]="3000 6000"
sizes[cas13]="2500 6000"

# These are the catalytic residues in the reference proteins
declare -A residues
residues[cas9]="D10 H840"
residues[cas12]="D908 E993"
residues[cas13]="R472 H477 R1048 H1053"



# =====================================
# Make sure the output directory exists
# =====================================

mkdir -p $OUTPUT
mkdir -p $REBLAST_OUTPUT_DIR



# ===========================================
# Find operons that meet our minimal criteria
# ===========================================

minimal_file=$OUTPUT/minimal_passing_operons.csv.gz
if [[ ! -e $minimal_file ]]; then
    for filename in $INPUT; do
        # stream tar files to stdout
        tar -xzOf $filename
    done | python fix_paths.py $KEEP_PATHS | python rules-all.py | python rules-nocas1-2.py | gzip > $minimal_file
else
    >&2 echo "Minimal file already exists. Skipping..."
fi 


# ======================================================================================================
# Apply our rules to the gene_finder output to find Class 1 systems. These have to be handled separately
# here and throughout the script because the way we determine whether they are nuclease dead is done 
# implicitly by the absence of Cas3 or Cas10.
# ======================================================================================================

class1_file=$OUTPUT/class1.csv.gz
if [[ ! -e $class1_file ]]; then
    gzip -cd $minimal_file | python rules-nocas3-10.py | python rules-class1.py | python dedup.py | gzip > $class1_file 
else
    >&2 echo "Class I operons already exist. Skipping..."
fi



# ============================================================================================
# For Class 2 systems, apply type-specific rules, deduplicate operons, and select systems with 
# appropriately-sized effectors
# ============================================================================================

for group in "${!sizes[@]}"; do
    if [[ ! -e "$OUTPUT/$group.csv.gz" ]]; then 
        protein_name="${proteins[$group]}"
        protein_range="${sizes[$group]}"
        echo $group $protein_name $protein_range
    else
        >&2 echo "$group operons already exist. Skipping..."
    fi
done | parallel --colsep ' ' "gzip -cd $minimal_file | python rules-{2}.py | python dedup.py | python size_select.py {2} {3} {4} | gzip > $OUTPUT/{1}.csv.gz"



# ============================================================================================
# Make a FASTA file of each effector protein. For Class I proteins, we're using Cas7 at it is 
# essential and does not appear as a fusion with Cas5/6/8. For Class 2 proteins, we'll just
# use the nuclease, either Cas9, Cas12 or Cas13. If more than one effector protein is present
# in an operon, we pick the one with the lowest e-value. We also make a separate file for the
# Class 2 systems, which has a reference protein at the top of the file. We will perform a
# multiple sequence alignment on this file to determine which proteins are nuclease-active.
# ============================================================================================

# Class 2 
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
# Class 1
if [[ ! -e $class1_clustering_file ]]; then
    gzip -cd $class1_file | python make_effector_fasta.py cas7 > $class1_clustering_file
fi



# ==================================================================================================
# Perform a multiple sequence alignment so we can identify nucleases with mutated catalytic residues
# ==================================================================================================

for group in "${!proteins[@]}"; do
    fasta=$OUTPUT/$group.for_alignment.fasta
    alignment=$OUTPUT/$group.afa
    mafft_strategy="$OUTPUT/$group.mafft_strategy"
    protein="${proteins[$group]}"
    residues_file="$OUTPUT/$group.residues.csv"
    nuclease_dead_operons="$OUTPUT/$group.nuclease_dead.csv.gz"

    # Perform the alignment
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



# =======================================================================================================
# For the nuclease-inactive systems (or Class 1 systems without a nuclease) we now group them together
# by protein similarity, as many contigs are highly similar to several others and this reduces the burden
# during the manual analysis, and also makes it easier to see patterns or deviations from patterns. 
# As before, we use Cas7 as the protein for Class 1 systems, and Cas9, Cas12 or Cas13 for Type 1 systems.
# =======================================================================================================

# Do the clustering
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

# If all of the operons in a cluster are nuclease-dead, we copy their data to a new directory for further processing
# We don't need to set a threshold as we saw empirically that there were zero clusters with both nuclease-active and
# nuclease dead systems
for group in "${!proteins[@]}"; do
    ndcluster_dir=$OUTPUT/$group.nuclease_dead_clusters
    if [[ -e $ndcluster_dir ]]; then
        echo "nuclease-dead analysis done for $group"
        continue
    fi
    mkdir -p $ndcluster_dir
    nuclease_dead_operons="$OUTPUT/$group.nuclease_dead.csv.gz"
    cluster_gzip_dir=$OUTPUT/$group.clusters
    for nuclease_dead_cluster_filename in $(python find_nuclease_dead_clusters.py $nuclease_dead_operons $cluster_gzip_dir); do
        cp $cluster_gzip_dir/$nuclease_dead_cluster_filename $ndcluster_dir
    done
done




# ===================================================================
# Find self-targeting spacers and inverted repeats in Class 2 systems
# =================================================================== 
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



# ===================================================================
# Find self-targeting spacers and inverted repeats in Class 1 systems
# =================================================================== 

fully_analyzed_dir=$OUTPUT/class1.fully-analyzed
if [[ ! -e $fully_analyzed_dir ]]; then
    mkdir -p $fully_analyzed_dir
    for filename in $(ls $OUTPUT/class1.clusters | rg '.*csv.gz'); do
        echo $filename
    done | parallel -j 8 --colsep ' ' "gzip -cd $OUTPUT/class1.clusters/{1} | python self_targeting.py | python find_inverted_repeats.py | gzip > $fully_analyzed_dir/{1}"
else
    >&2 echo "Inverted repeats and self-targeting spacers already found for Class 1"
fi



# =========================================================================
# Re-BLAST candidate systems with the TrEMBL, Swissprot, and tRNA databases
# =========================================================================

for group in "${!proteins[@]}"; do
    filedir=$OUTPUT/$group.fully-analyzed
    for filename in $(ls $filedir); do
        cluster_directory=$(basename $filename)
        directory="$REBLAST_OUTPUT_DIR/$group/$cluster_directory"
        mkdir -p $directory
        completed_file=$REBLAST_OUTPUT_DIR/$group/completed.txt

        if [[ ! $(rg $filename $completed_file) ]]; then
            >&2 echo "Reblasting $group operons from $filename"
            # Perform the re-BLASTing. If we skip the file because it has too few clusters, we print nothing so that we can lower the minimum cluster size later. If we do print the filename it indicates that we BLASTed as many operons as we wanted.
            # python reblast.py $BLASTN_DB $BLASTP_DB $MIN_REBLAST_CLUSTER_SIZE $REBLAST_COUNT $directory $OUTPUT/$group.clusters/$filename >> $completed_file
        else
            >&2 echo "Already finished $group operons from $filename"
        fi
    done
done



# =========================================================================
# Make side-by-side plots of the operons with annotations from our custom
# transposon/Cas databases (on top) and the annotations from the re-BLAST
# databases (on bottom)
# =========================================================================

for group in cas9 cas12 cas13 class1; do
    output_dir="$OUTPUT/plots/$group"
    if [[ ! -e $output_dir ]]; then
        >&2 echo "Plotting re-BLASTed operons for $group"
        fd '.*csv$' $OUTPUT/reblast/$group | parallel -j2 'cat {}' | python plot_reblast.py $OUTPUT/$group.fully-analyzed $output_dir
    else
        >&2 echo "Already plotted re-BLASTed operons for $group"
    fi
done



# =========================================================================================
# We manually looked at the reblasted results alongside their original
# annotations and decided which ones looked plausible enough to merit further investigation
# =========================================================================================
input_file=$OUTPUT/plausible-operon-ids.txt
output_file=$OUTPUT/plausible.csv.gz

if [[ -e $input_file && ! -e $output_file ]]; then
    fd '.*csv.gz' $OUTPUT/*fully-analyzed* | parallel -j2 'gzip -cd {}' | python extract_candidates.py $input_file | gzip > $output_file
fi
