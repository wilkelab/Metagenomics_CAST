""" Runs minCED on each system to find CRISPR arrays. For some reason, ~25% of Type V systems don't
get their arrays deteceted by piler-cr. """


import re
import subprocess
import sys
from uuid import uuid4
from operon_analyzer import genes, load

regex = re.compile(r"^(\d+)\s+(\w+)\s+(\w+).*")


def run_minced(operon):
    sequence = load.load_sequence(operon)
    opstart = max(0, operon.start-2000)
    opend = min(len(sequence), operon.end+2000)
    sequence = sequence[opstart: opend]
    uid = uuid4()

    with open(f"/tmp/mincedseq-{uid}.fa", "w") as f:
        f.write(f">seq\n{str(sequence)}")
    subprocess.call(f'minced -minNR 2 /tmp/mincedseq-{uid}.fa /tmp/minced-{uid}.crisprs'.split())
    with open(f"/tmp/minced-{uid}.crisprs") as f:
        for line in f:
            match = regex.match(line)
            if match:
                start = int(match.group(1))
                repeat = match.group(2)
                spacer = match.group(3)
                end = start + len(repeat) + len(spacer)
                feature = genes.Feature('Repeat Spacer', (start+opstart, end+opstart), '', None, '', None, f"{repeat} {spacer}", spacer)
                operon._features.append(feature)


if __name__ == '__main__':
    for operon in load.load_operons(sys.stdin):
        run_minced(operon)
        print(operon.as_str())
