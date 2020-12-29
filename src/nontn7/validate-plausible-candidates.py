# Creates a more readable text representation of each candidate system. The basic idea here is that
# we went through the operon diagrams and found things that seemed interesting, and then wanted to look up the 
# details to see what exactly a particular gene BLASTed to or whether CRISPR targets had reasonable homology
# and so forth.


import sys

from operon_analyzer import load, genes, visualize

candidate_csv_gz = sys.argv[1]


def print_details(operon: genes.Operon):
    for feature in operon:
        print(f"{feature.name}  {feature.start}..{feature.end}  {len(feature)} {feature.aln_len} {feature.pident}  {feature.qcovhsp}")


plausible_candidates = list(load.load_gzipped_operons(candidate_csv_gz))
reblasted_operons = list(load.load_operons(sys.stdin))


for plausible, reblasted in visualize.make_operon_pairs(plausible_candidates, reblasted_operons):
    print(f"{plausible.contig},{plausible.start}..{plausible.end}")
    print("========")
    print_details(reblasted)
    print("")
    for feature in plausible:
        if feature.name == 'CRISPR target':
            print(feature.name, feature.description)
        if feature.name.startswith('IR'):
            print(feature.name, feature.start, feature.end, feature.description)
        if feature.name == 'CRISPR array':
            print(feature.start, feature.end, feature.description, feature.sequence)
