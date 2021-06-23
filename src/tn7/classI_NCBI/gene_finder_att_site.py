from gene_finder.pipeline import Pipeline
import yaml
import os, sys

import os

in_path,out_path,tns_DB,cas_DB,att_DB,ffs_DB = sys.argv[1:]

def get_all_file_paths(path_to_dir):
    """
    Takes a path to a single directory and returns a list of the
    paths to each file in the directory. 
    """
    fasta_files = os.listdir(path_to_dir)
    full_paths = [os.path.join(path_to_dir, f) for f in fasta_files]
    return full_paths


def GF(name):
    """
    Command line args:
        1. Path to input fasta file
        2. Path to config file

    Constructs a Pipeline object to search for CRISPR-Transposon
    systems in the genome/config of interest.
    """
    min_prot_len= 80 # ORF aa cutoff length
    span=20000 # length (bp) upstream & downstream of bait hits to keep
    outfrmt="CSV" # output file format (can be either "CSV" or "JSON")
    gzip=False # make this True if the input contig file is gzipped
    genomic_data = name
    print("Blasting file {}.".format(name))
    
    
    # create a directory to write results to, if it doesn't exist already
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # results file will be called "<contig filename>_results.csv"
    job_id = os.path.basename(genomic_data)
    #outfile = os.path.join(conf["out-dir"], "{}_results.csv".format(job_id))

    p = Pipeline()

    p.add_seed_step(db=tns_DB, name='tns', e_val=0.001, blast_type="PROT")
    p.add_blastn_step(db=ffs_DB, name='ffs', e_val = 0.01, parse_descriptions=False)
    p.add_blast_step(db=cas_DB, name='cas', e_val=0.005, blast_type="PROT")
    p.add_blast_step(db=att_DB, name='att', e_val=0.001, blast_type="PROT")
    p.add_crispr_step()
    p.run(genomic_data, job_id, out_path, min_prot_len,span, gzip)


import pandas as pd
from multiprocessing import Pool
if __name__ == '__main__':

    all_fasta_files = get_all_file_paths(in_path)

    p = Pool(80)
    p.map(GF, all_fasta_files)
