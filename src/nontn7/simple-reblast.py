import sys

from operon_analyzer import load
from gene_finder import pipeline
from genbank import load_sequence


NUM_THREADS = 96


def find_system_bounds(operon, sequence):
    features = list(operon.all_features)
    positions = []
    for feature in features:
        positions.append(feature.start)
        positions.append(feature.end)

    left_buffer_end = min(positions)
    left_buffer_start = max(0, left_buffer_end - 3000)

    right_buffer_start = max(positions)
    right_buffer_end = min(len(sequence), right_buffer_start + 3000)
    return left_buffer_start, right_buffer_end


if __name__ == '__main__':
    blastn_db_dir = sys.argv[1]
    blastp_db_dir = sys.argv[2]
    output_dir = sys.argv[3]

    for n, operon in enumerate(load.load_operons(sys.stdin)):
        try:
            sequence = load_sequence(operon)
        except:
            sequence = load.load_sequence(operon)
        start, end = find_system_bounds(operon, sequence)
        job_id = f"{operon.contig}-{n}-{operon.start}-{operon.end}"
        p = pipeline.Pipeline()
        p.add_seed_with_coordinates_step(start=start,
                                         end=end,
                                         contig_id=operon.contig,
                                         db=blastp_db_dir,
                                         name="blastdb",
                                         e_val=1e-30,
                                         num_threads=NUM_THREADS,
                                         blast_type="PROT",
                                         parse_descriptions=False,
                                         blast_path='blastp-2.10')
        p.add_blastn_step(blastn_db_dir, 'tRNA', 1e-4, parse_descriptions=False, num_threads=NUM_THREADS, blastn_path='blastn-2.10')
        # run the pipeline on this one contig
        print(f"[RUNNING] Re-BLAST {job_id}", file=sys.stderr)
        results = p.run(data=operon.contig_filename, output_directory=output_dir, job_id=job_id, gzip=True)
