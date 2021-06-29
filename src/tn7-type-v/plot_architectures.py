""" Makes plots of each Type V CAST system, binning each architecture in separate directories. """

import os
import sys
from itertools import zip_longest

from operon_analyzer import load, visualize

architectures = {}
transposon_genes = set(['tnsA', 'tnsB', 'tnsC', 'tnsD', 'tnsE', 'tniQ', 'cas12k', 'ST', 'CANONSTS', 'STS', 'STR', 'PAM'])


plot_directory = sys.argv[1]


def plot_operons(operons, output_dir):
    for operon in operons:
        operon._features = [feature for feature in operon if feature.name not in ('STS source', 'STS target')]
        for feature in operon:
            if feature.name == 'Repeat Spacer':
                feature.description = ''
        try:
            filename = f"{operon.contig}-{operon.start}-{operon.end}.png"
            out_filename = os.path.join(output_dir, filename)
            fig = visualize.create_operon_figure(operon, False, color_by_blast_statistic='e_val', nucl_per_line=20000,
                                                 show_accession=False, show_description=True,
                                                 colormin=1e-30, colormax=1e-3)
            if fig is None:
                continue
            visualize.save_operon_figure(fig, out_filename)
        except Exception as e:
            print(e)
            continue


for operon in load.load_operons(sys.stdin):
    # First, find if any self-targeting spacers are canonical
    arrays = [feature for feature in operon if feature.name in ('Repeat Spacer', 'CRISPR array')]
    stss = operon.get('STS')
    for sts in stss:
        for array in arrays:
            if array.start < sts.start < array.end or array.start < sts.end < array.end:
                sts.name = 'CANONSTS'
                break

    features = []
    sorted_features = list(sorted(operon, key=lambda feat: feat.start))
    if len(sorted_features) < 2:
        continue
    deduplicated_features = []
    for feature1, feature2 in zip_longest(sorted_features, sorted_features[1:], fillvalue=None):
        if feature2 is not None and feature1.name == feature2.name and feature1.start == feature2.start and feature1.end == feature2.end and feature1.sequence == feature2.sequence:
            continue
        deduplicated_features.append(feature1)

    for feature in deduplicated_features:
        expected_gene = feature.name in transposon_genes
        good_annotation = feature.e_val is None or feature.e_val <= 1e-3
        if expected_gene and good_annotation:
            features.append(feature.name)

    features = tuple(features)
    if features in architectures:
        architectures[features].append(operon)
    elif tuple(reversed(features)) in architectures:
        architectures[tuple(reversed(features))].append(operon)
    else:
        architectures[features] = [operon]

# Now make plots, with each different architecture in its own directory.
for features, operons in architectures.items():
    subdirectory_name = '-'.join(features)
    directory = os.path.join(plot_directory, subdirectory_name)
    os.makedirs(directory, exist_ok=True)
    plot_operons(operons, directory)
