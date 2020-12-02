import gzip
import os
import random
import subprocess
import sys

from gene_finder import pipeline
from operon_analyzer import analyze

NUM_THREADS = 64

blastn_db_dir = sys.argv[1]
blastp_db_dir = sys.argv[2]
min_cluster_size = int(sys.argv[3])
reblast_count = int(sys.argv[4])
output_dir = sys.argv[5]
input_file = sys.argv[6]

with gzip.open(input_file, 'rt') as f:
    operons = sorted(analyze.load_operons(f), key=lambda operon: -len(operon))
if len(operons) < min_cluster_size:
    # This cluster has too few members, so we assume it's abberant and not worth looking into. This assumption could be
    # wrong! It's probably worth revisiting when the larger clusters have finished reblasting
    exit(0)

# Get the operon with the most features first
operons_to_reblast = [operons[0]]

# Get the operon with the most coverage, if it's not the same one
longest_operon = sorted([operon for operon in operons], key=lambda o: -abs(o.end-o.start))[0]
if longest_operon != operons_to_reblast[0]:
    operons_to_reblast.append(longest_operon)

remaining_operons = [op for op in operons if op not in operons_to_reblast]
if remaining_operons and reblast_count > len(operons_to_reblast):
    additional_operons = random.sample(remaining_operons, min(len(remaining_operons), reblast_count - len(operons_to_reblast)))
    operons_to_reblast.extend(additional_operons)

successful = 0
for operon in operons_to_reblast:
    assert os.path.exists(operon.contig_filename), f"missing file: {operon.contig_filename}"

    # set up the pipeline
    p = pipeline.Pipeline()
    p.add_seed_with_coordinates_step(start=operon.start,
                                     end=operon.end,
                                     contig_id=operon.contig,
                                     db=blastp_db_dir,
                                     name="blastdb",
                                     e_val=1e-30,
                                     num_threads=NUM_THREADS,
                                     blast_type="PROT",
                                     parse_descriptions=False,
                                     blast_path='blastp-2.10')
    p.add_crispr_step()
    p.add_blastn_step(blastn_db_dir, 'tRNA', 1e-30, parse_descriptions=False, num_threads=NUM_THREADS, blastn_path='blastn-2.10')
    # run the pipeline on this one contig
    try:
        job_id = f"{operon.contig}-{operon.start}-{operon.end}"
        results = p.run(data=operon.contig_filename, output_directory=output_dir, job_id=job_id, gzip=True)
    except subprocess.CalledProcessError:
        # pilercr segfaults for some unknown reason, in the event that this
        # happens we'll just skip this operon
        continue
    else:
        successful += 1
        if successful == reblast_count:
            break

# signal that we're done with this file
print(input_file)
