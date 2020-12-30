""" Takes the all_seqs results from mmseqs2 and a stream of Operons and saves each to a separate gzipped CSV """

import gzip
import os
import sys
from collections import defaultdict
from typing import IO, Dict, List, Tuple

from operon_analyzer import load

OperonLookup = Tuple[str, int, int]

cluster_file = sys.argv[1]
output_dir = sys.argv[2]


def _parse_cluster_file(text: IO) -> List[List[OperonLookup]]:
    # parses the mmseqs2 output and makes a list of lists of the members of each cluster.
    # the members are stored as the operon accession and start/end coordinates
    clusters = []

    for line in text:
        line = line.strip()
        if line.startswith(">"):
            nextline = next(text)
            if nextline.startswith(">"):
                # we're in a new group. line is the group, nextline is a member
                clusters.append([_parse_member(nextline[1:])])
            else:
                clusters[-1].append(_parse_member(line[1:]))
    return clusters


def _parse_member(member: str):
    # parse the unique sequence of each operon lookup
    accession, start, end, _, _, _ = member.split(",")
    start = int(start)
    end = int(end)
    return accession, start, end


def _make_sensible_groups(clusters: List[List[OperonLookup]]) -> Dict[OperonLookup, str]:
    # reformat the clusters so that we can check group membership easily while going through
    # the stream operon objects
    grouped = {}
    for n, cluster in enumerate(clusters):
        name = f'cluster{n}'
        for lookup in cluster:
            grouped[lookup] = name
    return grouped


if __name__ == '__main__':
    with open(cluster_file) as f:
        clusters = _parse_cluster_file(f)
        grouped = _make_sensible_groups(clusters)
        operons = load.load_operons(sys.stdin)
        operon_groups = defaultdict(list)
        for operon in operons:
            group = grouped.get((operon.contig, operon.start, operon.end))
            if not group:
                continue
            operon_groups[group].append(operon)
        for group, members in operon_groups.items():
            with gzip.open(os.path.join(output_dir, f'{group}.csv.gz'), 'w') as f:
                for member in members:
                    text = f"{member.as_str()}\n"
                    f.write(text.encode())
