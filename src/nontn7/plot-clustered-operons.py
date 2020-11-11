from operon_analyzer import analyze, visualize
from tools.filters import fs
import flabpal
import sys

image_dir = sys.argv[1]

feature_colors = {'transposase': flabpal.blue,
                  'cas1': flabpal.green,
                  'cas2': flabpal.blue,
                  'cas4': flabpal.blue,
                  'cas3': flabpal.red,
                  'cas9': flabpal.red,
                  'cas10': flabpal.red,
                  'cas12': flabpal.red,
                  'cas5': flabpal.pink,
                  'cas6': flabpal.brown,
                  'cas7': flabpal.purple,
                  'cas8': flabpal.yellow,
                  'CRISPR array': flabpal.orange,
                  '': flabpal.gray}


operons = list(analyze.load_operons(sys.stdin))
for operon in operons:
    fs.evaluate(operon)
clustered_operons = analyze.cluster_operons_by_feature_order(operons)

visualize._plot_clustered_operons(clustered_operons, image_dir, False, None, feature_colors=feature_colors, nucl_per_line=20000)
