""" Selects systems with CRISPR arrays located near Cas12k. """


import sys

from operon_analyzer import analyze, rules

from tools.filters import fs

rs = rules.RuleSet().max_distance('cas12k', 'Repeat Spacer', 1000, closest_pair_only=True)

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
