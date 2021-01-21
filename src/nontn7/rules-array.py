""" Accepts a stream of the raw gene_finder results and applies the most
basic possible rules. This should then be used to create a smaller dataset
that is more efficient to work with. """


import sys

from operon_analyzer import analyze, rules

from tools.filters import fs

rs = rules.RuleSet().require('CRISPR array')
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
