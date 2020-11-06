""" Searches for self-targeting spacers in operons, add them as Feature objects, and rewrite all operons back to stdout
regardless of what is found.
"""

from operon_analyzer import analyze, spacers
import sys


operons = analyze.load_operons(sys.stdin)
updated_operons = spacers.find_self_targeting_spacers(operons, 0.8, num_processes=32)
for operon in updated_operons:
    print(operon.as_str())
