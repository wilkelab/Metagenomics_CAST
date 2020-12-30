""" Computes values that we reference in the manuscript """

import re
import sys
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
from operon_analyzer import load

from tools.colors import blue, gray, orange, red

matplotlib.style.use('flab')


regexes = {'transposase': re.compile(r'transposase|transposition|transposon|integrase|integration|resolvase|recombinase|recombination|\bIS\d+|(T|t)np', re.IGNORECASE),
           'cas9': re.compile(r'\bcas9\w?', re.IGNORECASE),
           'cas12': re.compile(r'(\bcas12\w?)|(cpf1)', re.IGNORECASE),
           'cas13': re.compile(r'(\bcas13\w?)|(c2c2)', re.IGNORECASE),
           'cas6': re.compile(r'(\bcas6)|(\bcase)|(\bcsy4)', re.IGNORECASE),
           'cas5': re.compile(r'(\bcas5)|(\bcasd)|(\bcsc1)|(\bcsy2)|(\bcsf3)|(\bcsm4)|(\bcsx10)|(\bcmr3)', re.IGNORECASE),
           'cas8': re.compile(r'(\bcas8)|(\bcasa)|(\bcsh1)|(\bcsd1)|(\bcse1)|(\bcsy1)|(\bcsf1)', re.IGNORECASE),
           'cas7': re.compile(r'(\bcas7)|(\bcasc)|(\bcsd2)|(\bcsc2)|(\bcsy3)|(\bcsf2)|(\bcsm3)|(\bcsm5)|(\bcmr1)|(\bcmr6)|(\bcmr4)', re.IGNORECASE),
           'cas11': re.compile(r'(\bcas11)|(\bcasb)|(\bcse2)|(\bcsm2)|(\bcmr5)', re.IGNORECASE),
           'IR': re.compile(r'IR #\d+', re.IGNORECASE),
           'target': re.compile(r'CRISPR target'),
           'array': re.compile(r'(CRISPR )?array', re.IGNORECASE)}

counts = defaultdict(int)
facts = {'Cas only': 0, 'Tnp only': 0, 'Both': 0, 'Neither': 0}

for n, operon in enumerate(load.load_operons(sys.stdin)):
    operon_counts = defaultdict(int)
    for name, regex in regexes.items():
        for feature in operon:
            if regex.search(feature.name):
                operon_counts[name] += 1
    has_cas = False
    has_tnp = False
    has_ir = False
    has_target = False
    for protein in 'cas5', 'cas6', 'cas7', 'cas8', 'cas11':
        if operon_counts[protein] > 0:
            has_cas = True
            break
    has_tnp = operon_counts['transposase'] > 0
    has_ir = operon_counts['IR'] > 0
    has_target = operon_counts['target'] > 0

    if has_cas and has_tnp:
        facts['Both'] += 1
    elif has_cas and not has_tnp:
        facts['Cas only'] += 1
    elif not has_cas and has_tnp:
        facts['Tnp only'] += 1
    else:
        facts['Neither'] += 1

    for name, count in operon_counts.items():
        if count == 0:
            continue
        counts[name] += 1

print(facts)
fig, ax = plt.subplots()
labels = []
values = []
for name, count in sorted(facts.items(), key=lambda x: -x[1]):
    labels.append(name)
    values.append(count)

ax.bar([1, 2, 3, 4], values, color=[blue, red, orange, gray])
ax.set_xticklabels(labels)
ax.set_xticks([1, 2, 3, 4])
ax.set_ylabel("Number of putative operons")

plt.savefig("reannotation-results.png", bbox_inches='tight')


# * We found ZZZ systems with IRs alone and ZZZ systems with both IRs and TSDs.
# * We examined the spacer sequences of the systems and attempted to locate a corresponding target in the contig. We found ZZZ systems with targets of at least 80% homology that were within ZZZ base pairs of the IR or TSD.
# # * Current and historical gene names for cas (1-13) and Tn7 proteins (tnsA-E and tniQ, also glmS) were used as queries against the UniRef50 protein cluster databases (ZZZ URL?).
# * The three reference sets (tnsAB + transposases, tnsCDE + tniQ + glmS, and cas1-13) were each converted into blast databases using the NCBI blast+ command makeblastdb (version ZZZ) with default parameters. 
# * Approximately ZZZ% of the metagenomic operons that encoded a transposase and one Cas gene were nearly identical at the nucleotide sequence level. However, exact nucleotide comparisons were too slow to de-duplicate this large dataset. Instead, we considered two operons to be identical if they met the following properties: (1) they had the same protein-coding genes and CRISPR arrays in the same order; (2) the genes had the same relative distances to each other; and (3) the translated sequences of all proteins were identical. This de-duplication was applied to all operons prior to all downstream analysis.
# * We used MAFFT (Katoh and Standley, 2013) to align the Cas9, Cas12 and Cas13 protein sequences to each other, respectively. [zzz-Consider including a sentence on benchmarking this strategy to known Cas12k systems]. We found ZZZ% of Cas9, ZZZ% of Cas12, and ZZZ% of Cas13 had nuclease-inactivating mutations or deletions.


