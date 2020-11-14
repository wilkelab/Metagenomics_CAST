""" Determines the residues at a given position in an aligned sequence. The first sequence is a reference protein whose
catalytic residues are known. We then see what the equivalent residues are in all other proteins. This is later used to
determine which proteins might be nuclease-dead. 

"""

from Bio import SeqIO
import sys

# space-separated catalytic residues with their (1-based) index.
# For example, AsCas12a would be: D908 E993
catalytic_residues = sys.argv[1:]

# take user-supplied indices, converting to zero-indexed numbers
amino_acids = [aa[0] for aa in catalytic_residues]
interesting_indices = [int(val[1:]) - 1 for val in catalytic_residues]

# Load the MSA and get the reference protein's alignment
alignments = SeqIO.parse(sys.stdin, 'fasta')
reference = next(alignments)

# First, find the column of the catalytic residues in the reference protein
index = 0
indices = []
residues = []
for residue in reference:
    if residue != "-":
        residues.append(residue)
        indices.append(index)
    index += 1
columns = [indices[index] for index in interesting_indices]

# Make sure our user-supplied residues actually exist in the reference protein
for column, aa in zip(columns, amino_acids):
    assert reference[column] == aa

# Now write what each protein has in the catalytic protein position
for record in alignments:
    record_residues = "".join([record[column] for column in columns])
    print(f"{record_residues}\t{record.description}")
