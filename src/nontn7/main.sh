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
ARRAYS_OPTIONAL_DIRECTORY=$OUTPUT/systems-with-or-without-crispr-arrays
ARRAYS_REQUIRED_DIRECTORY=$OUTPUT/systems-with-crispr-arrays

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


# ======================================
# Make sure the output directories exist
# ======================================

for directory in $OUTPUT $REBLAST_OUTPUT_DIR $ARRAYS_REQUIRED_DIRECTORY $ARRAYS_OPTIONAL_DIRECTORY; do
    mkdir -p $directory
done

arrays_optional_file=$OUTPUT/minimal_passing_systems.csv.gz
arrays_required_file=$OUTPUT/minimal_passing_systems_with_crispr_arrays.csv.gz


function echodone () {
    >&2 echo "[DONE] $@"
}


function echodo () {
    >&2 echo "[RUNNING] $@"
}


function find_minimal_systems () {
    # Applies our bare-minimum rules to all of the neighborhoods in the gene_finder output
    # This produces two files: one where CRISPR arrays are optional, and one where they are required
    local action="Find minimal systems"

    if [[ ! -e $arrays_optional_file ]]; then
        echodo $action
        for filename in $INPUT; do
            # stream tar files to stdout
            tar -xzOf $filename
        done | python fix_paths.py $KEEP_PATHS | python rules-minimal.py | python rules-nocas1-2.py | gzip > $arrays_optional_file
    else
        echodone $action
    fi 

    # Make a separate dataset with the additional constraint that CRISPR arrays are required
    local action="Find minimal systems with CRISPR arrays"
    if [[ ! -e $arrays_required_file ]]; then
        echodo $action
        gzip -cd $arrays_optional_file | python rules-array.py | gzip > $arrays_required_file
    else
        echodone $action
    fi 
}


function find_minimal_subsystems () {
    # We select a subset of systems with some common trait to help focus our search.

    local systemname=$1
    local output="$OUTPUT/minimal-$systemname.csv.gz"
    local output_with_arrays="$OUTPUT/minimal-$systemname-with-arrays.csv.gz"

    local action="Find minimal systems with potential $systemname transposons"
    if [[ ! -e $output ]]; then
        echodo $action
        gzip -cd $arrays_optional_file | python rules-$systemname.py | python dedup.py | gzip > $output
    else
        echodone $action
    fi

    local action="Find minimal systems with potential $systemname transposons and CRISPR arrays"
    if [[ ! -e $output_with_arrays ]]; then
        echodo $action
        gzip -cd $output | python rules-array.py | gzip > $output_with_arrays
    else
        echodone $action
    fi
}


function apply_type_specific_rules () {
    # For Class 2 systems, apply type-specific rules, deduplicate operons, and select systems with 
    # appropriately-sized effectors

    local minimal_file=$1  # the raw input that contains all plausible systems
    local output_dir=$2  # directory to put results in
    mkdir -p $output_dir

    for group in "${!sizes[@]}"; do
        local action="Find $group systems for $output_dir from $minimal_file"
        if [[ ! -e "$output_dir/$group.csv.gz" ]]; then 
            echodo $action
            protein_name="${proteins[$group]}"
            protein_range="${sizes[$group]}"
            echo $group $protein_name $protein_range
        else
            echodone $action
        fi
    done | parallel --colsep ' ' "gzip -cd $minimal_file | python rules-{2}.py | python dedup.py | python size_select.py {2} {3} {4} | gzip > $output_dir/{1}.csv.gz"

    # ======================================================================================================
    # Apply our rules to the gene_finder output to find Class 1 systems. These have to be handled separately
    # here and throughout the script because the way we determine whether they are nuclease dead is done 
    # implicitly by the absence of Cas3 or Cas10.
    # ======================================================================================================

    class1_file=$output_dir/class1.csv.gz
    local action="Find Class 1 systems for $output_dir"
    if [[ ! -e $class1_file ]]; then
        echodo $action
        gzip -cd $minimal_file | python rules-nocas3-10.py | python rules-class1.py | python dedup.py | gzip > $class1_file 
    else
        echodone $action
    fi
}

function make_clustering_file () {
    # ============================================================================================
    # Make a FASTA file of each effector protein. For Class I proteins, we're using Cas7 at it is 
    # essential and does not appear as a fusion with Cas5/6/8. For Class 2 proteins, we'll just
    # use the nuclease, either Cas9, Cas12 or Cas13. If more than one effector protein is present
    # in an operon, we pick the one with the lowest e-value. We also make a separate file for the
    # Class 2 systems, which has a reference protein at the top of the file. We will perform a
    # multiple sequence alignment on this file to determine which proteins are nuclease-active.
    # ============================================================================================
    local directory=$1
    # Class 2 
    for group in "${!proteins[@]}"; do 
        protein="${proteins[$group]}"
        clustering_file="$directory/$group.for_clustering.fasta"
        local action="Make clustering file for $group"
        if [[ ! -e $clustering_file ]]; then
            echodo $action
            # Make the FASTA file to be used for clustering
            gzip -cd $directory/$group.csv.gz | python make_effector_fasta.py $protein > $clustering_file
            # Make the FASTA file to be used for alignment
            cat $DATA/nontn7/$group.fasta $clustering_file > $directory/$group.for_alignment.fasta
        else
            echodone $action
        fi
    done

    # Class 1
    class1_clustering_file="$directory/class1.for_clustering.fasta"
    if [[ ! -e $class1_clustering_file ]]; then
        echodo "Make clustering file for Class 1"
        gzip -cd $class1_file | python make_effector_fasta.py cas7 > $class1_clustering_file
    else
        echodone "Make clustering file for Class1"
    fi
}

perform_multiple_sequence_alignment () {
    # ==================================================================================================
    # Perform a multiple sequence alignment so we can identify nucleases with mutated catalytic residues
    # ==================================================================================================

    local directory=$1

    for group in "${!proteins[@]}"; do
        fasta=$directory/$group.for_alignment.fasta
        alignment=$directory/$group.afa
        mafft_strategy="$directory/$group.mafft_strategy"
        protein="${proteins[$group]}"
        residues_file="$directory/$group.residues.csv"
        nuclease_dead_operons="$directory/$group.nuclease_dead.csv.gz"

        # Perform the alignment
        local action="Perform MSA for $group"
        if [[ ! -e $alignment ]]; then
            echodo $action
            mafft --auto --thread 8 $fasta 2>$alignment.stderr > $alignment
        else
            echodone $action
        fi

        # Figure out what the strategy used by MAFFT was
        local action="Determine MAFFT strategy for $group"
        if [[ ! -e $mafft_strategy ]]; then
            strategy=$(rg --multiline 'Strategy:\n\s(.*)\n' $alignment.stderr --replace '$1') 
            echodo $action
            echo "$group $strategy" > $mafft_strategy
            rm $alignment.stderr
        else
            echodone $action
        fi

        # Use the alignment and determine which nucleases have (potentially) been inactivated
        local action="Identify catalytic residues for $group"
        if [[ ! -e $residues_file ]]; then
            echodo $action
            python identify_catalytic_residues.py ${residues[$protein]} < $alignment > $residues_file 
        else
            echodone $action
        fi

        # Gets operons with nuclease-dead effectors and save them to a new file
        local action="Separate systems with nuclease-dead effectors for $group"
        if [[ ! -e $nuclease_dead_operons ]]; then
            echodo $action
            gzip -cd $directory/$group.csv.gz | python load_nuclease_dead.py $protein $residues_file | gzip > $nuclease_dead_operons
        else
            echodone $action
        fi
    done
}

function cluster_nuclease_inactive_systems () {
    # For the nuclease-inactive systems (or Class 1 systems without a nuclease) we now group them together
    # by protein similarity, as many contigs are highly similar to several others and this reduces the burden
    # during the manual analysis, and also makes it easier to see patterns or deviations from patterns. 
    # As before, we use Cas7 as the protein for Class 1 systems, and Cas9, Cas12 or Cas13 for Type 1 systems.

    local directory=$1
    # Do the clustering
    for group in cas9 cas12 cas13 class1; do
        all_seqs_file="$directory/$group.clustered_all_seqs.fasta"
        clustering_file="$directory/$group.for_clustering.fasta"
        local action="Cluster $group systems by nuclease sequence"
        if [[ ! -e $all_seqs_file ]]; then
            echodo $action
            tempfile=/tmp/$group.mmseqs.temp
            mmseqs easy-cluster --min-seq-id 0.95 --write-lookup 1 $clustering_file $directory/$group.clustered $tempfile 2>/dev/null
            rm -r $tempfile
            mkdir -p $directory/$group.clusters
            gzip -cd $directory/$group.csv.gz | python get_clustered_operons.py $all_seqs_file $directory/$group.clusters
        else
            echodone $action
        fi
    done

    # If all of the operons in a cluster are nuclease-dead, we copy their data to a new directory for further processing
    # We don't need to set a threshold as we saw empirically that there were zero clusters with both nuclease-active and
    # nuclease dead systems
    for group in "${!proteins[@]}"; do
        ndcluster_dir=$directory/$group.nuclease_dead_clusters
        local action="Nuclease-dead analysis for $group"
        if [[ -e $ndcluster_dir ]]; then
            echodone $action
            continue
        fi
        echodo $action
        mkdir -p $ndcluster_dir
        nuclease_dead_operons="$directory/$group.nuclease_dead.csv.gz"
        cluster_gzip_dir=$directory/$group.clusters
        for nuclease_dead_cluster_filename in $(python find_nuclease_dead_clusters.py $nuclease_dead_operons $cluster_gzip_dir); do
            cp $cluster_gzip_dir/$nuclease_dead_cluster_filename $ndcluster_dir
        done
    done
}

function find_self_targeting_spacers_and_inverted_repeats () {

    local directory=$1
    # Find self-targeting spacers and inverted repeats in Class 2 systems
    for group in "${!proteins[@]}"; do
        fully_analyzed_dir=$directory/$group.fully-analyzed
        local action="Find inverted repeats and self-targeting spacers for $group"
        if [[ -e $fully_analyzed_dir ]]; then
            echodone $action 
            continue
        fi
        echodo $action
        mkdir -p $fully_analyzed_dir
        for filename in $(ls $directory/$group.nuclease_dead_clusters | rg '.*csv.gz'); do
            echo $group $filename
        done
    done | parallel -j 8 --colsep ' ' "gzip -cd $directory/{1}.nuclease_dead_clusters/{2} | python self_targeting.py | python find_inverted_repeats.py | gzip > $directory/{1}.fully-analyzed/{2}"

    # Find self-targeting spacers and inverted repeats in Class 1 systems
    fully_analyzed_dir=$directory/class1.fully-analyzed
    local action="Find inverted repeats and self-targeting spacers for Class 1"
    if [[ ! -e $fully_analyzed_dir ]]; then
        echodo $action
        mkdir -p $fully_analyzed_dir
        for filename in $(ls $directory/class1.clusters | rg '.*csv.gz'); do
            echo $filename
        done | parallel -j 8 --colsep ' ' "gzip -cd $directory/class1.clusters/{1} | python self_targeting.py | python find_inverted_repeats.py | gzip > $fully_analyzed_dir/{1}"
    else
        echodone $action
    fi
}

function reblast () {
    # Re-BLAST candidate systems with the TrEMBL, Swissprot, and tRNA databases
    local input_directory=$1
    for group in "${!proteins[@]}"; do
        filedir=$input_directory/$group.fully-analyzed
        for filename in $(ls $filedir); do
            python reblast.py $BLASTN_DB $BLASTP_DB $MIN_REBLAST_CLUSTER_SIZE $REBLAST_COUNT $REBLAST_OUTPUT_DIR $filedir/$filename
        done
    done
}

function plot_operons () {
    # Make side-by-side plots of the operons with annotations from our custom
    # transposon/Cas databases (on top) and the annotations from the re-BLAST
    # databases (on bottom)
    local input_directory=$1
    for group in cas9 cas12 cas13 class1; do
        output_dir="$input_directory/plots/$group"
        mkdir -p $output_dir
        >&2 echo "Plotting operons for $group in $input_directory"
        python plot_operons.py $input_directory/$group.fully-analyzed $output_dir
    done
}

function plot_reblasted_and_original_systems () {
    # Make side-by-side plots of the operons with annotations from our custom
    # transposon/Cas databases (on top) and the annotations from the re-BLAST
    # databases (on bottom)
    local input_directory=$1
    for group in cas9 cas12 cas13 class1; do
        output_dir="$input_directory/reblasted-plots/$group"
        >&2 echo "Plotting re-BLASTed operons for $group in $input_directory"
        if [[ ! -e $output_dir ]]; then
            mkdir -p $output_dir
            >&2 echo "Plotting re-BLASTed operons for $group in $input_directory"
            fd '.*csv$' $OUTPUT/reblast | parallel -j2 'cat {}' | python plot_reblast.py $input_directory/$group.fully-analyzed $output_dir
        else
            >&2 echo "Already plotted re-BLASTed operons for $group"
        fi
    done
}

function examine_plausible_candidates () {
    # We manually looked at the reblasted results alongside their original
    # annotations and decided which ones looked plausible enough to merit further investigation

    local input_file=$OUTPUT/plausible-operon-ids.txt
    local output_file=$OUTPUT/plausible.csv.gz

    if [[ -e $input_file && ! -e $output_file ]]; then
        fd '.*csv.gz' $OUTPUT/*fully-analyzed* | parallel -j2 'gzip -cd {}' | python extract_candidates.py $input_file | gzip > $output_file
    fi
}

function reblast_interesting_candidates () {
    local input_directory=$1
    for group in cas12 cas9 cas13 class1; do
        filedir=$input_directory/$group.fully-analyzed
        python reblast_interesting_candidates.py $BLASTN_DB $BLASTP_DB $REBLAST_COUNT $OUTPUT/interesting-candidates.csv $group $filedir $REBLAST_OUTPUT_DIR
    done
}    

function run_complete_analysis () {
    # Takes a single input file of some set of candidate systems, and runs the complete analysis for Cas9, Cas12, Cas13 and Class 1 systems.
    local input_systems=$1
    local output_directory=$2
    apply_type_specific_rules $input_systems $output_directory
    make_clustering_file $output_directory
    perform_multiple_sequence_alignment $output_directory
    cluster_nuclease_inactive_systems $output_directory
    find_self_targeting_spacers_and_inverted_repeats $output_directory
    plot_operons $output_directory
}

# Run the pipeline

find_minimal_systems
find_minimal_subsystems tn3
find_minimal_subsystems composite

# run_complete_analysis $arrays_required_file $ARRAYS_REQUIRED_DIRECTORY/all
# run_complete_analysis $OUTPUT/minimal-tn3-with-arrays.csv.gz $ARRAYS_REQUIRED_DIRECTORY/tn3
# run_complete_analysis $OUTPUT/minimal-composite-with-arrays.csv.gz $ARRAYS_REQUIRED_DIRECTORY/composite

# run_complete_analysis $OUTPUT/minimal-tn3.csv.gz $ARRAYS_OPTIONAL_DIRECTORY/tn3
# run_complete_analysis $OUTPUT/minimal-composite.csv.gz $ARRAYS_OPTIONAL_DIRECTORY/composite
# run_complete_analysis $arrays_optional_file $ARRAYS_OPTIONAL_DIRECTORY/all

reblast_interesting_candidates $ARRAYS_OPTIONAL_DIRECTORY/all
plot_reblasted_and_original_systems $ARRAYS_OPTIONAL_DIRECTORY/all
