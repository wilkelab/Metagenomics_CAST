# This deletes files in the output directory so that re-running the analysis can occur.
# It does NOT delete the minimal_passing_operons.csv.gz file or the first set of cas9/12/13/class1 csv.gz files.

OUTPUT=/stor/home/rybarskj/Metagenomics_CRISPR_transposons/output/nontn7

for key in mafft residues for_alignment cluster fully-analyzed afa nuclease_dead; do
    rm -rf $OUTPUT/*$key*
done

rm -rf $OUTPUT/plots/*
