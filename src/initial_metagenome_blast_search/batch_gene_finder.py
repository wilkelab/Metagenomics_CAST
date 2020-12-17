from gene_finder.pipeline import Pipeline
import yaml
import json, os, sys


genomic_data = sys.argv[1].strip() 
conf_file = sys.argv[2]
batch_name = sys.argv[3]
task_id = sys.argv[4]

# load contents of yaml conf file
stream = open(conf_file, 'r')
conf = yaml.safe_load(stream)
job_out = conf["out-dir"]
num_threads = conf["num-threads"]

# setup output directory, if it doesn't exist yet
if not os.path.exists(job_out):
    os.mkdir(job_out)

p = Pipeline()
for step in conf["steps"]:
        
    if step["type"] == "seed":
        p.add_seed_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"], num_threads=num_threads)
        
    elif step["type"] == "filter":
        p.add_filter_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"], num_threads=num_threads)
        
    elif step["type"] == "blast":
        p.add_blast_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"], num_threads=num_threads)
        
    else:
        p.add_crispr_step()

job_id = f"{batch_name}_{task_id}"
results = p.run(job_id=job_id, data=genomic_data, output_directory=job_out, min_prot_len=conf["min-prot-len"], span=conf["span"], gzip=conf["gzip"],
                record_all_hits=False, incremental_output=False, starting_contig=None)
