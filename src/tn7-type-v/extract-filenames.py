""" Gets the filenames for each system. This is just needed to make it easy to
parallelize the re-BLASTing step. """


import sys

from operon_analyzer import load


for operon in load.load_operons(sys.stdin):
    print(operon.contig_filename)
