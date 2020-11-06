import re
import pickle
import sys

from operon_analyzer import reannotation

# TODO: this needs to be totally rewritten. Use gzipped input instead of pickle files, and clean up the rest of the code after refactoring the summarization module in operon_analyze
passing_operons_pickle = sys.argv[1]
reblasted_operons_pickle = sys.argv[2]


with open(passing_operons_pickle, "rb") as f:
    raw_operons = pickle.load(f)

with open(reblasted_operons_pickle, "rb") as f:
    reblasted_operons = pickle.load(f)

osbad = re.compile(r"OS=.+")
otherbad = re.compile(r"\(.+\)")

for operon in reblasted_operons:
    for n, feature in enumerate(operon):
        feature.name = " ".join(feature.name.split()[1:])
        feature.name = osbad.sub("", feature.name)
        feature.name = otherbad.sub("", feature.name)

clusters, reblasted_operons = reannotation.prepare_operons_for_counting(raw_operons, reblasted_operons)
for label, cloperons in clusters.items():
    counts = reannotation.count_cluster_reannotations(cloperons, reblasted_operons)
    fractions = reannotation.convert_reannotation_counts_to_fractions(counts)

    if not counts:
        continue
    for (feature, count), reannotations in fractions.items():
        if feature != 'cas9':
            continue
        reannotated_labels = " ".join(reannotations.keys()).lower()
        if 'cas9' in reannotated_labels or 'nuclease' in reannotated_labels or 'CRISPR' in reannotated_labels:
            print(f"{len(cloperons)}-{'-'.join(label)}")
            # print(reannotation.format_reannotation_summary("-".join(label), len(cloperons), fractions))
            break
    else:
        continue
