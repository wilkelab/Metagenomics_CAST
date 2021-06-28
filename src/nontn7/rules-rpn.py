""" Selects systems with Cas12 and Rpn proteins. """


import sys

from operon_analyzer import analyze, rules

from tools.filters import fs

rs = rules.RuleSet().require('cas12', regex=True) \
                    .require('^Rpn', regex=True)

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
