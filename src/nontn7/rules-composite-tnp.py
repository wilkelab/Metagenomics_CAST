""" Finds systems that could potentially be composite transposons. """

import itertools
import sys

from operon_analyzer import analyze, genes, rules

from tools.filters import fs


def _might_be_composite(operon: genes.Operon) -> bool:
    # We consider systems with two transposases that BLAST to the same protein and have opposite orientation,
    # and at least one Cas gene is positioned between them, to be potentially composite
    transposases = [feature for feature in operon if feature.name == 'transposase']
    if len(transposases) < 2:
        return False
    cas_positions = [feature.start for feature in operon if feature.name.lower().startswith('cas')]
    if not cas_positions:
        return False

    pairs = itertools.combinations(transposases, 2)
    for left, right in pairs:
        if left.accession != right.accession or left.strand != right.strand:
            continue
        for cas_start in cas_positions:
            cas_is_contained_forward = left.start < cas_start < right.start
            cas_is_contained_reverse = left.start > cas_start > right.start
            if cas_is_contained_forward or cas_is_contained_reverse:
                return True
    return False


composite_rule = rules.Rule('composite', _might_be_composite)

rs = rules.RuleSet().custom(composite_rule).require('cas([5-9]|12|13)', regex=True)
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
