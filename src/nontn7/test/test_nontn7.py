import pytest

from make_effector_fasta import get_lowest_evalue_protein
from operon_analyzer.genes import Operon, Feature


def test_get_lowest_evalue_protein():
    feature1 = Feature('cas9', (3, 100), '', 1, '', 3e-90, '', 'MAAA', 100)
    feature2 = Feature('cas9', (200, 500), '', 1, '', 1e-30, '', 'MBBB', 100)
    operon = Operon('contig', '', 0, 10000, [feature1, feature2])

    protein = get_lowest_evalue_protein(operon, 'cas9')
    assert protein.sequence == 'MAAA'


def test_get_lowest_evalue_protein_missing():
    feature1 = Feature('cas9', (3, 100), '', 1, '', 3e-90, '', 'MAAA', 100)
    feature2 = Feature('cas9', (200, 500), '', 1, '', 1e-30, '', 'MBBB', 100)
    operon = Operon('contig', '', 0, 10000, [feature1, feature2])

    protein = get_lowest_evalue_protein(operon, 'cas12')
    assert protein is None


def test_get_lowest_evalue_protein_only_one():
    feature1 = Feature('transposase', (0, 100), '', 1, '', 1e-30, '', 'MAAA', 100)
    feature2 = Feature('cas9', (200, 500), '', 1, '', 1e-30, '', 'MBBB', 100)
    operon = Operon('contig', '', 0, 10000, [feature1, feature2])

    protein = get_lowest_evalue_protein(operon, 'cas9')
    assert protein.sequence == 'MBBB'
