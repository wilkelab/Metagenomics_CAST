import os
import random
import subprocess
import sys
from collections import defaultdict

from gene_finder import pipeline
from operon_analyzer import load

random.seed(43)

blastn_db_dir = sys.argv[1]
blastp_db_dir = sys.argv[2]
reblast_count = int(sys.argv[3])
interesting_candidate_filename = sys.argv[4]
protein_of_interest = sys.argv[5]
input_directory = sys.argv[6]
output_directory = sys.argv[7]

NUM_THREADS = 64


cluster_filenames = defaultdict(list)
with open(interesting_candidate_filename) as f:
    for line in f:
        protein, cluster_number = line.strip().split()
        cluster_filenames[protein].append(f"cluster{cluster_number}.csv.gz")

for filename in os.listdir(input_directory):
    if filename not in cluster_filenames[protein_of_interest]:
        continue
    operons = tuple(load.load_gzipped_operons(os.path.join(input_directory, filename)))
    # good_operons = []

    # Only look at systems with inverted repeats, since all the interesting systems had them
    # for operon in operons:
    #     if not any([name.startswith("IR #") for name in operon.feature_names]):
    #         continue
    #     good_operons.append(operon)

    good_operons = random.sample(operons, min(len(operons), reblast_count))
    for operon in good_operons:
        job_id = f"{operon.contig}-{operon.start}-{operon.end}"
        output_file = os.path.join(output_directory, f"{job_id}_results.csv")
        if os.path.exists(output_file):
            print(f"[DONE] Re-BLAST {protein_of_interest} {filename} {job_id}", file=sys.stderr)
        else:
            print(f"[RUNNING] Re-BLAST {protein_of_interest} {filename} {job_id}", file=sys.stderr)
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
                results = p.run(data=operon.contig_filename, output_directory=output_directory, job_id=job_id, gzip=True)
            except subprocess.CalledProcessError:
                # pilercr segfaults for some unknown reason, in the event that this
                # happens we'll just skip this operon
                continue
