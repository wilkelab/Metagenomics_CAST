import random
import logging
import os
import subprocess
import sys

from gene_finder import pipeline
from operon_analyzer import analyze

MIN_CLUSTER_SIZE = 5
REBLASTS_PER_CLUSTER = 1

# caution: main.sh processes all Class 2 datasets in parallel, actual thread usage will be 5x this number
NUM_THREADS = 16

blastn_db_dir = sys.argv[1]
blastp_db_dir = sys.argv[2]
output_dir = sys.argv[3]

protein = output_dir.split("/")[-1]

logging.basicConfig(filename=f"/tmp/reblast-{protein}.log", filemode='a', format='%(asctime)s     %(message)s', level=logging.DEBUG)
logging.info(f"protein group: {protein}")

operons = analyze.load_operons(sys.stdin)
clustered_operons = analyze.cluster_operons_by_feature_order(operons)
logging.info(f"{len(clustered_operons)} clusters destined for {output_dir}")

for label, cloperons in clustered_operons.items():
    if len(cloperons) < MIN_CLUSTER_SIZE:
        continue
    assert 'transposase' in label
    keydir = f"{len(cloperons)}-{'-'.join(label)}".replace("CRISPR array", "array")
    if len(keydir) > 200:
        # probably junk, skipping it
        continue
    cluster_output_dir = os.path.join(output_dir, keydir)
    if not os.path.exists(cluster_output_dir):
        os.makedirs(cluster_output_dir)
    if len(os.listdir(cluster_output_dir)) == REBLASTS_PER_CLUSTER:
        logging.info(f"{cluster_output_dir} is aleady done!")
        continue

    # Pick some operons at random, ensuring that they haven't been re-BLASTed already
    random.shuffle(cloperons)
    sample = []
    for operon in cloperons:
        job_id = f"{operon.contig}-{operon.start}-{operon.end}"
        output_file = f"{cluster_output_dir}/{job_id}.csv"
        if os.path.exists(output_file):
            continue
        else:
            sample.append(operon)
        if len(sample) == REBLASTS_PER_CLUSTER:
            break

    logging.info(f"Re-BLASTing {len(sample)} operons for {output_dir}/{keydir}")
    for operon in sample:
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
                                         parse_descriptions=False)
        p.add_crispr_step()
        p.add_blastn_step(blastn_db_dir, 'tRNA', 1e-30, parse_descriptions=False, num_threads=NUM_THREADS)

        job_id = f"{operon.contig}-{operon.start}-{operon.end}"
        assert os.path.exists(operon.contig_filename), f"missing file: {operon.contig_filename}"
        # run the pipeline on this one contig
        try:
            results = p.run(data=operon.contig_filename, output_directory=cluster_output_dir, job_id=job_id, gzip=True)
        except subprocess.CalledProcessError:
            # pilercr segfaults for some unknown reason, in the event that this
            # happens we'll just skip this operon
            continue
        else:
            logging.info(f"Finished BLASTING {job_id} for {cluster_output_dir}")
