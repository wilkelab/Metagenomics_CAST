import sys
from operon_analyzer import load, visualize
from tools.filters import fs
from tools.colors import feature_colors

output_dir = sys.argv[1]

operons = tuple(load.load_operons(sys.stdin))

for operon in operons:
    fs.evaluate(operon)
visualize.plot_operons(operons, output_dir, feature_colors=feature_colors, plot_ignored=False, nucl_per_line=20000)
