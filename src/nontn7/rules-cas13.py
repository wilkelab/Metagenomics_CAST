""" Selects all potential Type VI systems. """

import csv
import sys
from operon_analyzer import analyze, rules
from tools.filters import standard_fs as fs


csv.field_size_limit(sys.maxsize)
rs = rules.RuleSet().require('cas13')
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
