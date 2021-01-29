""" Makes operon diagrams with annotations from our curated database on top and TrEMBL/Swissprot/tRNA databases on bottom. """

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
cluster_regex = re.compile(r"cluster(\d+).csv.gz")

for operon in reblasted_operons:
    for n, feature in enumerate(operon):
        feature.name = " ".join(feature.name.split()[1:])
        feature.name = remove_os.sub("", feature.name)
        feature.name = remove_parentheticals.sub("", feature.name)

for original_operon_gz in os.listdir(original_operons_dir):
    cluster_number = int(cluster_regex.match(original_operon_gz).group(1))

    with gzip.open(os.path.join(original_operons_dir, original_operon_gz), 'rt') as f:
        original_operons = tuple(analyze.load_operons(f))
        for operon in original_operons:
            fs.evaluate(operon)

        assert original_operons, "Invalid original operons file."
        try:
            for operon, other in visualize.make_operon_pairs(original_operons, reblasted_operons):
                filename = f"cluster{cluster_number:03d}-{operon.contig}-{operon.start}-{operon.end}.png"
                out_filename = os.path.join(output_dir, filename)
                visualize.plot_operon_pair(operon, other, None, None, out_filename, None, False, feature_colors)
        except ValueError as e:
            print(e)
            print(f"Couldn't plot {os.path.join(original_operons_dir, original_operon_gz)}")
            continue
