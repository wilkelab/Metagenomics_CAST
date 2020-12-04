import gzip
import os
import re
import sys

from operon_analyzer import analyze, visualize

from tools.colors import feature_colors
from tools.filters import fs

original_operons_dir = sys.argv[1]
output_dir = sys.argv[2]

reblasted_operons = tuple(analyze.load_operons(sys.stdin))

# The Uniprot/Swissprot annotations have a bunch of information
# that makes it impossible to see the actual protein name when
# plotting operons, so we remove those.
remove_os = re.compile(r"OS=.+")
remove_parentheticals = re.compile(r"\(.+\)")

for operon in reblasted_operons:
    for n, feature in enumerate(operon):
        feature.name = " ".join(feature.name.split()[1:])
        feature.name = remove_os.sub("", feature.name)
        feature.name = remove_parentheticals.sub("", feature.name)

for original_operon_gz in os.listdir(original_operons_dir):
    print(original_operon_gz)
    cluster_output_dir = os.path.join(output_dir, original_operon_gz)
    os.makedirs(cluster_output_dir, exist_ok=True)

    with gzip.open(os.path.join(original_operons_dir, original_operon_gz), 'rt') as f:
        original_operons = tuple(analyze.load_operons(f))
        for operon in original_operons:
            fs.evaluate(operon)

        assert original_operons, "Invalid original operons file."
        try:
            visualize.plot_operon_pairs(original_operons, reblasted_operons, cluster_output_dir, feature_colors=feature_colors)
        except ValueError:
            print(f"Couldn't plot {os.path.join(original_operons_dir, original_operon_gz)}")
            continue
