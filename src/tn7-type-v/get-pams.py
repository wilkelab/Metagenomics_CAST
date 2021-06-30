# Writes the sequences for the direct repeat, target and PAM for each self-target in the Type V systems
# This produces a flat text file that will later be parsed to determine the sequence logo for the PAM
# in Figure 5. We do this so we can account for situations where the pairwise alignment was off by a
# base or more. This is a reasonable assumption when the repeat ends in GAAAG, which has both experimental
# and bioinformatic evidence to support it as the true end of the repeat. For other repeat sequences, we
# suspect they are misidentified self-targets, and in reality those systems do not target the genome.

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
