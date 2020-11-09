""" Selects potential Class I systems. """


import sys
from operon_analyzer import analyze, rules
from tools.filters import fs



rs = rules.RuleSet().contains_at_least_n_features(['cas5', 'cas6', 'cas7', 'cas8'], 3)
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
