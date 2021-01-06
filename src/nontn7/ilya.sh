set -euo pipefail 

OGOUT=../../output/nontn7
OUT=../../output/nontn7/ilya-composite-tn3


declare -A proteins
proteins[cas9]="cas9"
proteins[cas12]="cas12"
proteins[cas13]="cas13"

# Sizes are the minimum and maximum length of the effector genes that we allow, in base pairs
declare -A sizes
sizes[cas9]="2000 6000"
sizes[cas12]="3000 6000"
sizes[cas13]="2500 6000"


for group in "${!sizes[@]}"; do
    if [[ ! -e "$OUT/$group.csv.gz" ]]; then 
        protein_name="${proteins[$group]}"
        protein_range="${sizes[$group]}"
        >&2 echo "$group "
        echo $group $protein_name $protein_range
    else
        >&2 echo "$group operons already exist. Skipping..."
    fi
done | parallel --colsep ' ' "gzip -cd $OGOUT/composite.csv.gz | python rules-{2}.py | python size_select.py {2} {3} {4} | gzip > $OUT/{2}.csv.gz"

>&2 echo "class1"
gzip -cd $OGOUT/composite.csv.gz | python rules-class1.py | gzip > $OUT/class1.csv.gz

for group in cas9 cas12; do 
    clustering_file="$OUT/$group.for_clustering.fasta"
    if [[ ! -e $clustering_file ]]; then
        gzip -cd $OUT/$group.csv.gz | python make_effector_fasta.py $group > $clustering_file
    fi
done

clustering_file="$OUT/class1.for_clustering.fasta"
if [[ ! -e $clustering_file ]]; then
    gzip -cd $OUT/class1.csv.gz | python make_effector_fasta.py cas7 > $clustering_file
fi


for group in cas9 cas12 class1; do
    all_seqs_file="$OUT/$group.clustered_all_seqs.fasta"
    clustering_file="$OUT/$group.for_clustering.fasta"
    if [[ ! -e $all_seqs_file ]]; then
        tempfile=/tmp/$group.mmseqs.temp
        mmseqs easy-cluster --min-seq-id 0.95 --write-lookup 1 $clustering_file $OUT/$group.clustered $tempfile
        rm -r $tempfile
        mkdir -p $OUT/$group.clusters
        gzip -cd $OUT/$group.csv.gz | python get_clustered_operons.py $all_seqs_file $OUT/$group.clusters
    fi
done

