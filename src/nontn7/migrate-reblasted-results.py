# This is a "one-time use" script that moves re-BLAST results into the correct directories. The cluster IDs are
# randomized with each run so if I re-run my analysis, I need to migrate them to the new directory where they belong.
# Re-BLASTing takes an extremely long time relative to everything else so this is unavoidable.

import sys
import os
from operon_analyzer import load

# Where the gzipped CSV files from the new analysis:
# Probably output/nontn7/$protein.fully_analyzed
current_operon_dir = sys.argv[1]

# The directory to make new cluster directories and to place the reblasted operons in:
# Probably output/nontn7/reblast/$protein

output_dir = sys.argv[2]


reblasted_operons = {f'{operon.contig}-{operon.start}-{operon.end}_results.csv': operon
                     for operon in load.load_operons(sys.stdin)}

for gzipped_csv in os.listdir(current_operon_dir):
    if not gzipped_csv.endswith(".csv.gz"):
        continue
    current_path = os.path.join(current_operon_dir, gzipped_csv)
    current_operons = load.load_gzipped_operons(current_path)
    directory = os.path.join(output_dir, gzipped_csv)
    for operon in current_operons:
        filename = f'{operon.contig}-{operon.start}-{operon.end}_results.csv'
        reblasted = reblasted_operons.get(filename)
        if reblasted:
            os.makedirs(directory, exist_ok=True)
            output_filepath = os.path.join(directory, filename)
            with open(output_filepath, 'w') as f:
                f.write(reblasted.as_str())
