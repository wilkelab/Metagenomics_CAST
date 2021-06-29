""" Selects all potential Type V systems. """


import sys

from operon_analyzer import analyze, rules

import custom_rules
from tools.filters import fs

tn7_proteins = ("tnsA", "tnsB", "tnsC", "tnsD", "tniQ")

rs = rules.RuleSet().custom(rule=custom_rules.contains_at_least_n_features_n_bp_apart(feature_list=tn7_proteins,
                                                                                      feature_count=2,
                                                                                      distance_bp=1500,
                                                                                      count_multiple_copies=True))

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
