""" Searches for self-targeting spacers in operons, add them as Feature objects, and rewrite all operons back to stdout regardless of what is found. """

import sys

from operon_analyzer import analyze, spacers

from tools.filters import fs

operons = analyze.load_operons(sys.stdin)
filtered_operons = []
for operon in operons:
    fs.evaluate(operon)
    filtered_operons.append(operon)
updated_operons = spacers.find_self_targeting_spacers(filtered_operons, 0.8, num_processes=32)
for operon in updated_operons:
    print(operon.as_str())
