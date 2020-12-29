""" Scans the entire contig where each putative operon is located and searches for inverted repeats """


import sys

from operon_analyzer import analyze, repeat_finder

from tools.filters import fs

for operon in analyze.load_operons(sys.stdin):
    fs.evaluate(operon)
    repeat_finder.find_inverted_repeats(operon, 500, 15)
    print(operon.as_str())
