""" Counts the number of operons that are passed into the script from stdin and prints the result to stdout """

import sys

from operon_analyzer import load

print(sum([1 for operon in load.load_operons(sys.stdin)]))
