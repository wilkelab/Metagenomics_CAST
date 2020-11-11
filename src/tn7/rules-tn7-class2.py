""" Selects minimal class I systems with nearby Tn7 features """

import csv
import sys
from operon_analyzer import analyze, rules
import custom_rules

csv.field_size_limit(sys.maxsize)

fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9)
rs = rules.RuleSet().custom(rule=custom_rules.contains_at_least_n_features_n_bp_apart(feature_list=["tnsA", "tnsB", "tnsC", "tnsD", "tniQ"], feature_count=2, distance_bp=1500, count_multiple_copies=False)) \
                    .max_distance(feature1_name=r'cas(9|12\w?|13\w?|phi)', feature2_name="CRISPR array", distance_bp=10000, closest_pair_only=True, regex=True)

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
