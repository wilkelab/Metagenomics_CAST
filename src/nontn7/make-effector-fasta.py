# Creates a FASTA file containing the protein sequences for a given protein.
# This is used to create input for MAFFT and ultimately to identify nuclease-dead
# proteins.

from operon_analyzer import analyze
from tools.filters import fs
import sys

protein_name = sys.argv[1]

# get the reference protein that we'll use to identify the catalytic residues with
reference_protein_filename = f"fasta/{protein_name}.fasta"
with open(reference_protein_filename) as f:
    reference_text = f.read()

sequences = []
for operon in analyze.load_operons(sys.stdin):
    fs.evaluate(operon)
    for n, feature in enumerate(operon.get(protein_name)):
        sequences.append((operon.contig, operon.start, operon.end, feature.accession, feature.sequence, n))

print(reference_text)
for contig, start, end, accession, sequence, n in sequences:
    print(f">{contig} {start} {end} {accession} {n}")
    print(f"{sequence}")
