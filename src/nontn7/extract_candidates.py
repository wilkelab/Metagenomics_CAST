# Takes a list of candidate systems that was produced during
# the manual validation, finds the gene_finder output for those systems,
# and prints their serialized data to stdout


import sys
from operon_analyzer import load


candidate_list = sys.argv[1]


with open(candidate_list) as f:
    candidates = set()
    for line in f:
        accession, coords = line.strip().split(',')
        coords = coords.split("..")
        start, end = int(coords[0]), int(coords[1])
        candidates.add((accession, start, end))


for operon in load.load_operons(sys.stdin):
    if (operon.contig, operon.start, operon.end) in candidates:
        print(operon.as_str())
