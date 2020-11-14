import re
import sys

from operon_analyzer import reannotation, load


operons = load.load_gzipped_operons(sys.argv[1])
reblasted_operons = tuple(load.load_operons(sys.stdin))

# Some regular expressions to clean up some irrelevant text from Uniprot gene names
osbad = re.compile(r"OS=.+")
otherbad = re.compile(r"\(.+\)")

for operon in reblasted_operons:
    for n, feature in enumerate(operon):
        feature.name = " ".join(feature.name.split()[1:])
        feature.name = osbad.sub("", feature.name)
        feature.name = otherbad.sub("", feature.name)

reannotation.summarize(operons, reblasted_operons)


# clusters, reblasted_operons = reannotation.prepare_operons_for_counting(operons, reblasted_operons)
# for label, cloperons in clusters.items():
#     counts = reannotation.count_cluster_reannotations(cloperons, reblasted_operons)
#     fractions = reannotation.convert_reannotation_counts_to_fractions(counts)

#     if not counts:
#         continue
#     for (feature, count), reannotations in fractions.items():
#         if feature != 'cas9':
#             continue
#         reannotated_labels = " ".join(reannotations.keys()).lower()
#         if 'cas9' in reannotated_labels or 'nuclease' in reannotated_labels or 'CRISPR' in reannotated_labels:
#             print(f"{len(cloperons)}-{'-'.join(label)}")
#             # print(reannotation.format_reannotation_summary("-".join(label), len(cloperons), fractions))
#             break
#     else:
#         continue
