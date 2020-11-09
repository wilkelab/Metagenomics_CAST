""" Requires that systems have cas1 and cas2. """


import sys
from operon_analyzer import analyze, rules
from tools.filters import fs




rs = rules.RuleSet().require('cas1').require('cas2')

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
