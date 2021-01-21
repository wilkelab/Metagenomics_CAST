""" Finds systems that have Tn3-family transposases. """

import sys

from operon_analyzer import analyze, genes, rules

from tools.filters import fs


def _has_tn3(operon: genes.Operon) -> bool:
    for feature in operon:
        if 'tn3' in feature.description.lower():
            return True
    return False


tn3_rule = rules.Rule('composite', _has_tn3)

rs = rules.RuleSet().custom(tn3_rule)
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
