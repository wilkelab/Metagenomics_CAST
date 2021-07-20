""" Gets sequences for our patent application. """

from operon_analyzer import load
from operon_analyzer import rules

import os
from tools.filters import tn7fs

import sys

relevent_genes = ('tnsB', 'tnsC', 'tniQ', 'cas12k')
tn7_genes = ('tnsB', 'tnsC', 'tniQ')


def get_pattern_features(operon, pattern):
    """ Pulls out a set of Features that have certain genes in a certain order.
    Returns None if not found. Returns all matching sets. """
    tn7fs.evaluate(operon)
    features = sorted([feature for feature in operon if feature.name in relevent_genes], key=lambda feature: feature.start)
    if len(features) < len(pattern):
        return []

    if len(features) == len(pattern):
        return [list(check_for_match(features, pattern))]

    zippers = []
    matches = []
    for n, _ in enumerate(pattern):
        zippers.append(features[n:])

    for subset_features in zip(*zippers):
        match = check_for_match(subset_features, pattern)
        if match:
            matches.append(list(match))
    return matches


def check_for_match(features, pattern):
    assert len(features) == len(pattern)

    # make sure we're not just looking at genes that blasted to more than one thing
    coords = set([(feature.start, feature.end) for feature in features])
    if len(coords) < len(pattern):
        return []

    # now do the comparisons to see if we matched the given pattern
    forward = all([feature.name == pat for feature, pat in zip(features, pattern)])
    reverse = all([feature.name == pat for feature, pat in zip(features, reversed(pattern))])
    if forward:
        return list(features)
    elif reverse:
        return list(reversed(features))
    else:
        return []


def get_matches_with_close_tn7_genes(matches, max_distance=300):
    for features in matches:
        if not features:
            continue
        for f1, f2 in zip(features, features[1:]):
            if f1.name not in tn7_genes or f2.name not in tn7_genes:
                continue
            if rules._feature_distance(f1, f2) > max_distance:
                break
        else:
            yield features


def run(operon, pattern):
    sequence = load.load_sequence(operon)
    matches = get_pattern_features(operon, pattern)
    good_matches = list(get_matches_with_close_tn7_genes(matches))
    cas12s = [feature for feature in operon if feature.name == 'cas12k']

    for match_number, match in enumerate(good_matches):
        orientations = set([feature.strand for feature in match])
        if len(orientations) != 1:
            continue

        coords = []
        for m in match:
            coords.append(m.start)
            coords.append(m.end)
        start = min(coords)
        end = max(coords)

        orientation = list(orientations)[0]
        closest_cas12 = None
        best_distance = None
        if orientation == 1:
            for cas12 in cas12s:
                if cas12.start > end:
                    distance = cas12.start - end
                    if best_distance is None or distance < best_distance:
                        closest_cas12 = cas12
                        best_distance = distance
        elif orientation == -1:
            for cas12 in cas12s:
                if cas12.end < start:
                    distance = start - cas12.end
                    if best_distance is None or distance < best_distance:
                        closest_cas12 = cas12
                        best_distance = distance
        if closest_cas12 is None:
            continue
        match.append(closest_cas12)
        coords.append(closest_cas12.start)
        coords.append(closest_cas12.end)

        start = min(coords)
        end = max(coords)

        if orientation == 1:
            start = max(0, start - 200)
            end = min(len(sequence), end + 700)
        elif orientation == -1:
            start = max(0, start - 700)
            end = min(len(sequence), end + 200)

        nucleotide_sequence = sequence[start: end]
        names = [feature.name for feature in match]
        # fix the mislabelled tnsC. we need to go back and perform an MSA to make sure this is correct in all cases
        if 'tnsC' not in names:
            match[1].name = 'tnsC'

        protein_fasta = "\n".join([f">{feature.name}\n{feature.sequence}" for feature in match])
        nucleotide_fasta = f">{operon.contig}_{operon.start}_{operon.end}\n{nucleotide_sequence}"

        yield match_number, operon, protein_fasta, nucleotide_fasta


if __name__ == '__main__':
    data_dir = sys.argv[1]
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(os.path.join(data_dir, 'protein'), exist_ok=True)
    os.makedirs(os.path.join(data_dir, 'nucleotide'), exist_ok=True)

    for operon in load.load_operons(sys.stdin):
        found = False
        for pattern in (('tnsB', 'tnsC', 'tniQ'),
                        ('tnsB', 'tnsB', 'tniQ'),
                        ('tnsB', 'tnsB', 'tnsC', 'tniQ'),
                        ('tnsB', 'tnsC', 'tnsB', 'tnsC', 'tniQ')):
            for match_number, operon, protein_fasta, nucleotide_fasta in run(operon, pattern):
                found = True
                with open(f"{data_dir}/protein/{operon.contig}_{operon.start}_{operon.end}_{match_number}.protein.fa", "w") as f:
                    f.write(protein_fasta)
                with open(f"{data_dir}/nucleotide/{operon.contig}_{operon.start}_{operon.end}_{match_number}.fa", "w") as f:
                    f.write(nucleotide_fasta)
        if not found:
            print(operon.as_str())
