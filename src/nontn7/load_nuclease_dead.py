# Creates a gene_finder-formatted operon CSV, containing only operons with putatively dead nucleases.

from operon_analyzer import analyze
import re
import sys

protein = sys.argv[1]
residues_file = sys.argv[2]

nuclease_active_residue_regex = {"cas9": re.compile(r'(D|E)(R|H|K)'),
                                 "cas12": re.compile(r'(D|E){2}'),
                                 "cas13": re.compile(r'(R|H|K){4}')}
regex = nuclease_active_residue_regex[protein]

# note which operons still contain a catalytically-active CRISPR-Cas nuclease, so that we can exclude them later
# many operons have multiple hits to our protein of interest so we can't make a determination by just examining
# the first one.
nuclease_active = set()
with open(residues_file) as f:
    for line in f:
        motif, operon_data = line.strip().split()
        if regex.match(motif):
            contig_accession, contig_start, contig_end, _, _, _ = operon_data.split(',')
            nuclease_active.add((contig_accession, int(contig_start), int(contig_end)))

# write the nuclease-dead operons to stdout
for operon in analyze.load_operons(sys.stdin):
    if (operon.contig, operon.start, operon.end) in nuclease_active:
        continue
    sys.stdout.write(operon.as_str())
