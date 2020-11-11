import gzip
import re
import sys

from operon_analyzer import analyze, visualize
from tools.colors import feature_colors

original_operons_gz = sys.argv[1]
output_dir = sys.argv[2]

reblasted_operons = list(analyze.load_operons(sys.stdin))

osbad = re.compile(r"OS=.+")
otherbad = re.compile(r"\(.+\)")

for operon in reblasted_operons:
    for n, feature in enumerate(operon):
        feature.name = " ".join(feature.name.split()[1:])
        feature.name = osbad.sub("", feature.name)
        feature.name = otherbad.sub("", feature.name)

with gzip.open(original_operons_gz, 'rt') as f:
    original_operons = list(analyze.load_operons(f))

visualize.plot_operon_pairs(original_operons, reblasted_operons, output_dir, feature_colors=feature_colors)
