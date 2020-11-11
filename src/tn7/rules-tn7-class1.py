""" Selects candidates with a class I effector module and a group of at least 2 Tn7 genes """

import csv
import sys
from operon_analyzer import analyze, rules
import custom_rules

csv.field_size_limit(sys.maxsize)

fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9)
rs = rules.RuleSet().custom(rule=custom_rules.contains_at_least_n_features_n_bp_apart(feature_list=["tnsA", "tnsB", "tnsC", "tnsD", "tniQ"], feature_count=2, distance_bp=1500, count_multiple_copies=False)) \
                    .custom(rule=custom_rules.contains_at_least_n_features_n_bp_apart(feature_list=["cas5", "cas6", "cas7", "cas8"], feature_count=3, distance_bp=1500, count_multiple_copies=False)) 

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
