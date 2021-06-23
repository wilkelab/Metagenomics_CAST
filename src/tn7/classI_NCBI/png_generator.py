import sys
from operon_analyzer.analyze import load_analyzed_operons
from operon_analyzer.visualize import build_operon_dictionary, plot_operons

analysis_csv, pipeline_csv, image_directory = sys.argv[1:]
good_operons = []

with open(pipeline_csv) as f:
    operons = build_operon_dictionary(f)
with open(analysis_csv) as f:
    for contig, filename, start, end, result in load_analyzed_operons(f):
        if result[0] != 'pass':
            continue
        op = operons.get((contig, filename, start, end))
        if op is None:
            continue
        good_operons.append(op)
print("{} figures have been generated.".format(len(good_operons)))
for x in good_operons:
    try:
        plot_operons([x], image_directory)
        print('1')
    except:
        raise
