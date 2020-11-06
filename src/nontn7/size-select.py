""" Selects systems with a Feature that is within a given size range. """

import csv
import sys
from operon_analyzer import analyze, rules
from tools.filters import standard_fs as fs

feature_name = sys.argv[1]
min_bp = int(sys.argv[2])
max_bp = int(sys.argv[3])

assert min_bp <= max_bp, "Lower size limit must be less than or equal to the upper limit"

csv.field_size_limit(sys.maxsize)

rs = rules.RuleSet().minimum_size(feature_name, min_bp) \
                    .maximum_size(feature_name, max_bp)

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
