""" Counts how many systems came from which of our two databases. """

from operon_analyzer import load

import sys

ncbi = 0
embl = 0
other = 0

for operon in load.load_operons(sys.stdin):
    if 'NCBI' in operon.contig_filename:
        ncbi += 1
    elif 'EMBL_EBI' in operon.contig_filename:
        embl += 1
    else:
        other += 1

print(f"NCBI: {ncbi}")
print(f"EMBL: {embl}")
print(f"Other: {other}")
print(f"Total:{ncbi+embl+other}")
