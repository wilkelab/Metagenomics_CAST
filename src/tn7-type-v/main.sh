#!/usr/bin/env bash
# This script runs the Tn7 Type V CRISPR-transposon pipeline.


set -euo pipefail  # Halt the pipeline if any errors are encountered

OUTPUT=../../output/tn7-type-v
DATA=../../data/tn7-type-v
INPUT=$(ls /stor/work/Wilke/amh7958/pipeline-results/*tar.gz /stor/work/Wilke/amh7958/pipeline-results/missing_effectors/*.tar.gz)
DATA=../../data/tn7-type-v
KEEP_PATHS="NO"

mkdir -p $OUTPUT $OUTPUT/plots $OUTPUT/reblast

if [[ ! -e $OUTPUT/cas12k.csv.gz ]]; then
    >&2 echo "Running systems through Cas12k rules"
    for filename in $INPUT; do
        # stream tar files to stdout
        tar -xzOf $filename
    done | python fix_paths.py NO | python rules-cas12k.py | gzip > $OUTPUT/cas12k.csv.gz
fi

if [[ ! -e $OUTPUT/cas12k-tns.csv.gz ]]; then
    >&2 echo "Running systems through Tn7 rules"
    gzip -cd $OUTPUT/cas12k.csv.gz | python rules-tn7.py | python dedup.py | gzip > $OUTPUT/cas12k-tns.csv.gz
fi

if [[ ! -e $OUTPUT/cas12k-tns-minced.csv.gz ]]; then
    >&2 echo "Finding CRISPR arrays"
    gzip -cd $OUTPUT/cas12k-tns.csv.gz | python minced.py | gzip > $OUTPUT/cas12k-tns-minced.csv.gz
fi

if [[ ! -e $OUTPUT/cas12k-tns-minced-array.csv.gz ]]; then
    >&2 echo "Filtering systems without arrays"
    gzip -cd $OUTPUT/cas12k-tns-minced.csv.gz | python rules-array.py | gzip > $OUTPUT/cas12k-tns-minced-array.csv.gz 
fi

if [[ ! -e $OUTPUT/unique-filenames.csv ]]; then
    >&2 echo "Getting unique filenames"
    gzip -cd $OUTPUT/cas12k-tns-minced-array.csv.gz | python extract-filenames.py | sort -u > $OUTPUT/unique-filenames.csv
fi

mkdir -p $OUTPUT/reblast

if [[ "$(ls $OUTPUT/reblast | wc -l)" -eq "0" ]]; then
    >&2 echo "Re-BLASTing systems"
    cat $OUTPUT/unique-filenames.csv | parallel -j 96 "python relook.py {} $OUTPUT/reblast"
fi

mkdir -p $OUTPUT/reblast-minced
if [[ "$(ls $OUTPUT/reblast-minced | wc -l)" -eq "0" ]]; then
    >&2 echo "Running minced on re-BLASTed results"
    for filename in $(ls $OUTPUT/reblast/*csv); do
        basename $filename
    done | parallel -j 96 "cat $OUTPUT/reblast/{} | python minced.py > $OUTPUT/reblast-minced/{}"
fi

if [[ ! -e $OUTPUT/reblasted-cas12k-with-sts.csv.gz ]]; then
    >&2 echo "Finding self-targeting spacers"
    for filename in $(ls $OUTPUT/reblast-minced/*csv); do
        echo $filename
    done | parallel -j 96 "python sts.py {}" | gzip > $OUTPUT/reblasted-cas12k-with-sts.csv.gz
fi

if [[ ! -e $OUTPUT/reblasted-cas12k-with-sts-rules-deduped.csv.gz ]]; then
    gzip -cd $OUTPUT/reblasted-cas12k-with-sts.csv.gz | python rules-cas12k.py | python rules-tn7.py | python rules-array.py | python dedup.py | gzip > $OUTPUT/reblasted-cas12k-with-sts-rules-deduped.csv.gz
fi

# Do the work needed to determine the target site identities

if [[ ! -e $OUTPUT/self-target-context-seqs.fa ]]; then
    >&2 echo "Making fasta file of target sequences so we can BLAST them"
    gzip -cd  $OUTPUT/reblasted-cas12k-with-sts.csv.gz | python make-self-target-fasta.py > $OUTPUT/self-target-context-seqs.fa
fi

if [[ ! -e $OUTPUT/self-target-blastn-nt.csv ]]; then
    >&2 echo "blastn att sites"
    blastn-2.10 -db ~/work/nt/nt -query $OUTPUT/self-target-context-seqs.fa -evalue 1e-8 -outfmt '6 qseqid sseqid evalue qseq sseq pident stitle' -num_threads 16 -max_target_seqs 100 | rg -v 'genome' | rg -v 'chromosome' | rg -v 'complete sequence' | rg -v 'genomic( DNA)? sequence' > $OUTPUT/self-target-blastn-nt.csv
fi

if [[ ! -e $OUTPUT/self-target-blastx-nr.csv ]]; then
    >&2 echo "blasttx att sites"
    blastx-2.10 -db ~/work/nrdb/nr -query $OUTPUT/self-target-context-seqs.fa -evalue 1e-8 -outfmt '6 qseqid sseqid evalue qseq sseq pident stitle' -num_threads 16 -max_target_seqs 100 | rg -v 'genome' | rg -v 'chromosome' | rg -v 'complete sequence' | rg -v 'genomic( DNA)? sequence' > $OUTPUT/self-target-blastx-nr.csv
fi

# Perform alignments for Figure 5

if [[ ! -e $OUTPUT/fig5d.afa ]]; then
    >&2 echo "Fig 5D alignments"
    mafft --auto $DATA/fig-5d-tnsc.fa > $OUTPUT/fig5d.afa
fi

if [[ ! -e $OUTPUT/fig5e.afa ]]; then
    >&2 echo "Fig 5E alignments"
    mafft --auto $DATA/fig-5e-tnsb.fa > $OUTPUT/fig5e.afa
fi

if [[ ! -e $OUTPUT/canonical-pams.csv ]]; then
    >&2 echo "Getting sequences of self-targeting PAMs"
    gzip -cd $OUTPUT/reblasted-cas12k-with-sts-rules-deduped.csv.gz | python get-pams.py > $OUTPUT/canonical-pams.csv 2> $OUTPUT/other-pams.csv
fi
