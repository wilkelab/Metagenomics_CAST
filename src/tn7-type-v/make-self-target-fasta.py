""" Makes FASTA files of the DNA 50 bp around each target so we can BLAST them. """

import sys

from operon_analyzer import load


for operon in load.load_operons(sys.stdin):
    sequence = load.load_sequence(operon)
    for feature in operon:
        if feature.name == 'ST':
            start = max(0, feature.start - 50)
            end = min(len(sequence), feature.end + 50)
            target_sequence = sequence[start: end]
            print(f">{operon.contig}_{operon.start}..{operon.end}_{feature.start}..{feature.end}")
            print(str(target_sequence))
