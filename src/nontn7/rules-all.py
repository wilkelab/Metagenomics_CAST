""" Accepts a stream of the raw gene_finder results and applies the most
basic possible rules. This should then be used to create a smaller dataset
that is more efficient to work with. """


import sys

from operon_analyzer import analyze, rules

from tools.filters import fs

rs = rules.RuleSet().require('CRISPR array') \
                    .require('transposase') \
                    .require(r'cas(5|6|7|8|9|10|11|12\w?|13\w?|phi)', regex=True) \
                    .exclude('tns[A-E]', regex=True) \
                    .exclude('tniQ')

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
