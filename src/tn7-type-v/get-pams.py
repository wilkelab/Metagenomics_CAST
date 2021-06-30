# Writes the sequences for the PAM for each self-target in the Type V systems.
# This produces a flat text file that will be used to generate a sequence logo for the PAM
# in Figure 5. There are some cases where the pairwise alignment seems to have made an
# off-by-one error - in these cases, we adjust the alignment by one base pair. This is a
# reasonable assumption given that there is both bioinformatic and experimental evidence
# that the atypical repeat always ends in "GAAAG".
# We split the output into two groups - those systems with canonical atypical repeats,
# and all others. We believe the latter group is the result of spurious alignments from
# systems that don't actually have self-targeting spacers (or are using them in trans).


import sys

from operon_analyzer import load

for operon in load.load_operons(sys.stdin):
    data = {"STS": [],
            "STR": [],
            "PAM": [],
            "ST": []}
    for feature in operon:
        if feature.name in ('STS', 'STR', 'PAM', 'ST'):
            data[feature.name].append(feature)

    seen_pams = []
    for spacer, repeat, pam, target in zip(data["STS"], data["STR"], data["PAM"], data["ST"]):

        # If two systems target the exact same location, we don't want to double count that PAM
        if pam.start in seen_pams:
            continue
        seen_pams.append(pam.start)

        if repeat.sequence.endswith("GAAA") and spacer.sequence.startswith("G"):
            # We've found a case that's probably an off-by-one error, so we adjust the sequences
            repeat.sequence = repeat.sequence[1:] + "G"
            spacer.sequence = spacer.sequence[1:]
            pam.sequence = pam.sequence[1:] + target.sequence[0]
            target.sequence = target.sequence[1:]

        if repeat.sequence.endswith("GAAAG"):
            outfile = sys.stdout
        else:
            outfile = sys.stderr
        print(f"{pam.sequence}", file=outfile)
