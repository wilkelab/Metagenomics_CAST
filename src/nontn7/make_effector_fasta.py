""" Prints FASTA-formatted text containing the sequences of the given protein, from the operons passed into stdin.
If there is more than one protein of the same name in the same operon, we pick the one with the lowest evalue. """

import sys
from typing import Optional

from operon_analyzer import genes, load

from tools.filters import fs


def get_lowest_evalue_protein(operon: genes.Operon, protein_name: str) -> Optional[genes.Feature]:
    proteins = sorted([feature for feature in operon.get(protein_name, regex=True)], key=lambda x: x.e_val)
    if proteins:
        return proteins[0]
    return None


if __name__ == '__main__':
    protein_name = sys.argv[1]

    sequences = []
    for operon in load.load_operons(sys.stdin):
        fs.evaluate(operon)
        protein = get_lowest_evalue_protein(operon, protein_name)
        if protein:
            sequences.append((operon.contig, operon.start, operon.end, protein.accession, protein.start, protein.end, protein.sequence))

    for contig, start, end, accession, pstart, pend, sequence in sequences:
        print(f">{contig},{start},{end},{accession},{pstart},{pend}")
        print(f"{sequence}")
