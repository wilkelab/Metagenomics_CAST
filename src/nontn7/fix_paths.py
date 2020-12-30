""" In the raw gene_finder output, paths to the FASTA files are valid on TACC (where gene_finder was run), but do not
point to the correct location on our local cluster. This script updates the paths in the Operon objects so that we can
access raw nucleotide data in downstream analyses. """


import csv
import os
import sys

from operon_analyzer import parse


def fix_path(contig_filename: str) -> str:
    # update paths from those at TACC to those on our local cluster
    leaf = contig_filename.replace("/scratch/07227/hilla3/projects/CRISPR-Transposons/data/genomic/", "")
    basedir = '/stor/scratch/Wilke/contig/database/metagenomic_contig_database'
    if not leaf.startswith("NCBI"):
        basedir = f'{basedir}/EMBL_EBI'
    good_path = os.path.join(basedir, leaf)
    assert os.path.exists(good_path)
    return good_path


if __name__ == '__main__':
    skip = sys.argv[1]
    if skip == 'YES':
        # This is useful is the pipeline is being run on a system other than our own
        for line in sys.stdin:
            print(line)
    else:
        writer = csv.writer(sys.stdout, delimiter=',')
        for line in parse.read_pipeline_output(sys.stdin):
            line[21] = fix_path(line[21])
            writer.writerow(line)
