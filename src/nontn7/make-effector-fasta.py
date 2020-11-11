# Creates a FASTA file containing the protein sequences for a given protein.
# This is used to create input for MAFFT and ultimately to identify nuclease-dead
# proteins.

from operon_analyzer import analyze
from tools.filters import fs
import sys

protein_name = sys.argv[1]
# TODO: fix these hardcoded paths
if len(sys.argv) > 2:
    reference_protein_filename = f"../../data/nontn7/{sys.argv[2]}.fasta"
else:
    reference_protein_filename = f"../../data/nontn7/{protein_name}.fasta"

# get the reference protein that we'll use to identify the catalytic residues with
with open(reference_protein_filename) as f:
    reference_text = f.read()

sequences = []
for operon in analyze.load_operons(sys.stdin):
    fs.evaluate(operon)
    for n, feature in enumerate(operon.get(protein_name)):
        sequences.append((operon.contig, operon.start, operon.end, feature.accession, feature.sequence, n))

print(reference_text.strip())
for contig, start, end, accession, sequence, n in sequences:
    print(f">{contig} {start} {end} {accession} {n}")
    print(f"{sequence}")
