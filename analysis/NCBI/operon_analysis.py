import sys
from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet, FilterSet



rs = RuleSet().require('tnsB',regex=True).require('cas7',regex=True).require('cas6',regex=True).require('tnsA',regex=True).require('glmS').require('tniQ').max_distance('tnsA', 'tnsB', 2000,closest_pair_only=True).max_distance('glmS', 'tniQ', 2000,closest_pair_only=True)

fs = FilterSet().pick_overlapping_features_by_bit_score(0.9).must_be_within_n_bp_of_anything(20000)

with open('NCBI.csv', 'w', newline='') as csvfile:
    with open('NCBI_genefinder_merge.csv', "r") as f:
        analyze(f, rs, fs,csvfile)
