# Searches an input contig (in fasta format) for regions with at least one cas gene (1-13) and one transposase no more than 25 kbp apart 

from gene_finder.pipeline import Pipeline
import yaml
import json, os, sys


genomic_data = sys.argv[1] # path to input file, in fasta format
output_directory = sys.argv[2] # directory to write output to

# setup output directory, if it doesn't exist
if not os.path.exists(output_directory):
    os.mkdir(output_directory)

p = Pipeline()
p.add_seed_step(db="data/initial_metagenomic_blast_search/databases/transposase", name="transposase", e_val=0.001, blast_type="PROT", num_threads=1)
p.add_filter_step(db="data/initial_metagenomic_blast_search/databases/cas_all", name="cas_all", e_val=0.001, blast_type="PROT", num_threads=1)
p.add_blast_step(db="data/initial_metagenomic_blast_search/databases/tn7-accessory", name="tn7_accessory", e_val=0.001, blast_type="PROT", num_threads=1)
p.add_crispr_step()

# use the input filename as the job id
# results will be written to the file <job id>_results.csv
job_id = os.path.basename(genomic_data)
p.run(job_id=job_id, data=genomic_data, output_directory=output_directory, min_prot_len=90, span=25000)
