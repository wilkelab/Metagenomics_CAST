""" Re-runs the search for Type V CASTs, but uses the Cas12k database as the seed instead of the transposase one.
We do this because we noticed that some Type V CASTs contain co-localized transposons or were actually insertions
at the exact same attachment site. """


import os
import sys

from gene_finder.pipeline import Pipeline

filename = sys.argv[1]
out_path = sys.argv[2]
min_prot_len = 80  # ORF aa cutoff length
span = 100000  # length (bp) upstream & downstream of bait hits to keep


filename = filename.strip()
job_id = os.path.basename(filename)
p = Pipeline()
p.add_seed_step(db='/stor/work/Wilke/amh7958/databases/cas12k/blast_db', name='cas12k', e_val=1e-30, blast_type="PROT", blast_path='blastp', num_threads=1)
p.add_blast_step(db='/stor/work/Wilke/blastDB/classifer/tns', name='tns', e_val=1e-3, blast_type='PROT', blast_path='blastp', num_threads=1)
p.add_crispr_step()
p.run(filename, job_id, out_path, min_prot_len, span, False, gzip=True)
