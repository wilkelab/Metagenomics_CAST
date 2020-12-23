from gene_finder.pipeline import Pipeline
import yaml
import os, sys

import os
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
    
    genomic_data = name
    print(name)
    conf_file = "./config.yaml"

    # load contents of yaml conf file
    stream = open(conf_file, 'r')
    conf = yaml.load(stream)
    
    # create a directory to write results to, if it doesn't exist already
    if not os.path.exists(conf["out-dir"]):
        os.mkdir(conf["out-dir"])

    # results file will be called "<contig filename>_results.csv"
    job_id = os.path.basename(genomic_data)
    #outfile = os.path.join(conf["out-dir"], "{}_results.csv".format(job_id))

    p = Pipeline()
    for step in conf["steps"]:
        
        if step["type"] == "seed":
            p.add_seed_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
        
        elif step["type"] == "filter":
            p.add_filter_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
        
        elif step["type"] == "blast":
            p.add_blast_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
        
        else:
            p.add_crispr_step()
    
    p.run(data=genomic_data, job_id=job_id, output_directory=conf["out-dir"], min_prot_len=conf["min-prot-len"],
          span=conf["span"], gzip=conf["gzip"])


import pandas as pd
from multiprocessing import Pool
if __name__ == '__main__':

    all_fasta_files = get_all_file_paths('./manualfasta')
    print(len(all_fasta_files))

    #p = Pool(50)
    p = Pool(70)
    p.map(GF, all_fasta_files)
