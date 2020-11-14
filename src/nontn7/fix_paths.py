""" In the raw gene_finder output, paths to the FASTA files are valid on TACC (where gene_finder was run), but do not
point to the correct location on our local cluster. This script updates the paths in the Operon objects so that we can
access raw nucleotide data in downstream analyses.
"""

import csv
import sys

from operon_analyzer import parse
from tools.diskio import fix_path


if __name__ == '__main__':
    writer = csv.writer(sys.stdout, delimiter=',')
    for line in parse.read_pipeline_output(sys.stdin):
        line[21] = fix_path(line[21])
        writer.writerow(line)
