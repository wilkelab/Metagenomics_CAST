""" Removes systems with cas3 and cas10. """


import sys

from operon_analyzer import analyze, rules

from tools.filters import fs

rs = rules.RuleSet().exclude('cas3').exclude('cas10')
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
