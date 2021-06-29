""" Selects all potential Type V systems. """


import sys

from operon_analyzer import analyze, rules

from tools.filters import fs


rs = rules.RuleSet().require('cas12k')


analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
