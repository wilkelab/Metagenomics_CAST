""" Selects systems with Cas12 and Rpn proteins. """


import sys

from operon_analyzer import analyze, genes, rules

from tools.filters import fs


def _has_rpn_rule(operon: genes.Operon):
    for feature in operon:
        if feature.description.startswith('Rpn'):
            return True
    return False


rpn_rule = rules.Rule('has-rpn', _has_rpn_rule)

rs = rules.RuleSet().require('cas12', regex=True) \
                    .custom(rpn_rule)

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
