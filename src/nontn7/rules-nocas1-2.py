""" Excludes systems that have Cas1 or Cas2. """

import sys

from operon_analyzer import analyze, rules

from tools.filters import fs

rs = rules.RuleSet().exclude('cas1').exclude('cas2')
analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)
