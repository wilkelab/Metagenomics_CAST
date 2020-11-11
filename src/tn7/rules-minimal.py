""" 
Selects candidates that meet minimum criteria for the Tn7 analysis: 

    1. At least one Cas effector gene must be present
    2. At least one gene from the Tn7 core machinery must be present
    3. A TniQ family protein (either TniQ or TnsD, in our pipeline) must be present

Candidates that pass this set of rules will be used as the input 
for the class1 or class2- specific analysis in the next step of the pipeline. 
"""

import csv
import sys
from operon_analyzer import analyze, rules

csv.field_size_limit(sys.maxsize)

fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9)
rs = rules.RuleSet().require(r'cas(5|6|7|8|9|12\w?|13\w?|phi)', regex=True) \
                    .require(r'tns(A|B|C)', regex=True) \
                    .require("tniQ|tnsD", regex=True)

analyze.evaluate_rules_and_reserialize(sys.stdin, rs, fs)