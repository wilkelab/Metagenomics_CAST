""" Takes a list of candidate systems that was produced during
the manual validation, finds the gene_finder output for those systems,
and prints their serialized data to stdout """


import sys

from operon_analyzer import load

candidate_list = sys.argv[1]


with open(candidate_list) as f:
    candidates = []
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        _, accession, start, end = line.split('-')
        start, end = int(start), int(end)
        candidates.append((accession, start, end))


found = 0
operons = tuple(load.load_operons(sys.stdin))
for accession, start, end in candidates:
    for operon in operons:
        if (operon.contig, operon.start, operon.end) == (accession, start, end):
            found += 1
            print(operon.as_str())
            print("\n\n")
            break

if found != len(candidates):
    raise ValueError("Not all specified candidate systems were found!")
