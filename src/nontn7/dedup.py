# Filters out putative operons that are extremely similar to one another.

from operon_analyzer import analyze
import sys


operons = analyze.load_operons(sys.stdin)
unique = analyze.deduplicate_operons_approximate(operons)
unique = analyze.dedup_supersets(unique)

for operon in unique:
    print(operon.as_str())
