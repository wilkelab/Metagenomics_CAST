""" Selects all potential Type II systems. """

import sys
from operon_analyzer import analyze, rules
from tools.filters import fs


rs = rules.RuleSet().require('cas9')
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
