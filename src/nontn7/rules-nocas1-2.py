""" Removes systems with cas1 or cas2. """

import csv
import sys
from operon_analyzer import analyze, rules
from tools.filters import standard_fs as fs


csv.field_size_limit(sys.maxsize)
rs = rules.RuleSet().exclude('cas1').exclude('cas2')
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
