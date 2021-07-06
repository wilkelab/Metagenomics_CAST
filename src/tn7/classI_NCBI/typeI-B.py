import sys
from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet, FilterSet


input_genefinder_path,out_path = sys.argv[1:]

rs = RuleSet().require('tnsB').require('tnsA').require('cas7').require('cas6').require('cas5').require('glmS').max_distance('tnsA', 'tnsB', 1000,closest_pair_only=True).max_distance('tnsB', 'glmS', 3000,closest_pair_only=True)

fs = FilterSet().pick_overlapping_features_by_bit_score(0.9).must_be_within_n_bp_of_anything(10000)

with open(out_path, 'w', newline='') as csvfile:
    with open(input_genefinder_path, "r") as f:
        analyze(f, rs, fs,csvfile)
