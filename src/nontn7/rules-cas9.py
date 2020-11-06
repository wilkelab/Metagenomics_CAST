""" Selects all potential Type II systems. """

import csv
import sys
from operon_analyzer import analyze, rules
from tools.filters import standard_fs as fs, no_scarlet_mutation_filter


csv.field_size_limit(sys.maxsize)


rs = rules.RuleSet().require('cas9')
fs = fs.custom(no_scarlet_mutation_filter)
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
