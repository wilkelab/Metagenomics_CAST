""" Separates protein sequences for making phylogenetic trees. We also keep the low e-value
proteins because TnsC is usually mislabelled as TnsB. """

import sys

from operon_analyzer import load

protein = sys.argv[1]
limit = float(sys.argv[2])

evalues = []
good_seqs = set()
bad_seqs = set()

for operon in load.load_operons(sys.stdin):
    n = 0
    m = 0
    for feature in operon:
        if feature.name != protein:
            continue
        if feature.e_val < limit and feature.sequence not in good_seqs:
            print(f">{operon.contig}_{operon.start}_{operon.end}_{feature.start}_{feature.end}_{n}")
            print(f"{feature.sequence}")
            good_seqs.add(feature.sequence)
            n += 1
        elif feature.e_val >= limit and feature.sequence not in bad_seqs:
            print(f">{operon.contig}_{operon.start}_{operon.end}_{feature.start}_{feature.end}_{n}", file=sys.stderr)
            print(f"{feature.sequence}", file=sys.stderr)
            bad_seqs.add(feature.sequence)
            m += 1
