""" Helper functions for doing I/O. """

import os


def fix_path(contig_filename: str) -> str:
    # update paths from those at TACC to those on our local cluster
    leaf = contig_filename.replace("/scratch/07227/hilla3/projects/CRISPR-Transposons/data/genomic/", "")
    basedir = '/stor/scratch/Wilke/contig/database/metagenomic_contig_database'
    if not leaf.startswith("NCBI"):
        basedir = f'{basedir}/EMBL_EBI'
    good_path = os.path.join(basedir, leaf)
    assert os.path.exists(good_path)
    return good_path
