OUT=/stor/home/rybarskj/Metagenomics_CRISPR_transposons/output/nontn7/ilya-composite-tn3

for protein in cas9 cas12 class1; do
    cluster_gzs=$(ls $OUT/$protein.clusters);
    for filename in $cluster_gzs; do
        mkdir -p $OUT/plots/$protein/$filename
        gzip -cd $OUT/$protein.clusters/$filename | python plot_operons.py $OUT/plots/$protein/$filename
    done
done
