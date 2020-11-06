""" Selects potential Class I systems. """

import csv
import sys
from operon_analyzer import analyze, rules
from tools.filters import standard_fs as fs


csv.field_size_limit(sys.maxsize)
rs = rules.RuleSet().contains_at_least_n_features(['cas5', 'cas6', 'cas7', 'cas8'], 3)
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
