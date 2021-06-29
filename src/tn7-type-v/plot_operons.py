import os
import sys
from operon_analyzer import visualize, load

from tools.colors import feature_colors

output_dir = sys.argv[1]
os.makedirs(output_dir, exist_ok=True)


for operon in load.load_operons(sys.stdin):
    try:
        filename = f"{operon.contig}-{operon.start}-{operon.end}.png"
        out_filename = os.path.join(output_dir, filename)
        fig = visualize.create_operon_figure(operon, False, feature_colors, color_by_blast_statistic=False, nucl_per_line=20000, show_accession=False, show_description=True)
        if fig is None:
            continue
        visualize.save_operon_figure(fig, out_filename)
    except Exception as e:
        print(e)
        continue
