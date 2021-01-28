import gzip
import os
import re
import sys

from operon_analyzer import analyze, visualize

from tools.colors import feature_colors
from tools.filters import fs

cluster_regex = re.compile(r"cluster(\d+).csv.gz")


original_operons_dir = sys.argv[1]
output_dir = sys.argv[2]

for original_operon_gz in os.listdir(original_operons_dir):
    cluster_number = int(cluster_regex.match(original_operon_gz).group(1))
    with gzip.open(os.path.join(original_operons_dir, original_operon_gz), 'rt') as f:
        for operon in analyze.load_operons(f):
            fs.evaluate(operon)

            # don't plot systems without inverted repeats or self-targeting spacers
            has_ir_or_sts = any([feature for feature in operon if feature.name.startswith("IR #") or feature.name.startswith("CRISPR target")])
            if not has_ir_or_sts:
                continue
            try:
                filename = f"cluster{cluster_number:03d}-{operon.contig}-{operon.start}-{operon.end}.png"
                out_filename = os.path.join(output_dir, filename)
                fig = visualize.create_operon_figure(operon, False, feature_colors, color_by_blast_statistic=False, nucl_per_line=20000, show_accession=False, show_description=True)
                if fig is None:
                    continue
                visualize.save_operon_figure(fig, out_filename)
            except:
                continue
