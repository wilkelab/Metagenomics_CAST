""" Removes systems with cas1 or cas2. """


import sys
from operon_analyzer import analyze, rules
from tools.filters import fs



rs = rules.RuleSet().exclude('cas1').exclude('cas2')
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
