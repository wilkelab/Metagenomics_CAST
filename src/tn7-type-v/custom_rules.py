""" Extra rules needed to select Type V CASTs. """


from operon_analyzer.rules import Rule, RuleSet, _feature_distance
from operon_analyzer.genes import Operon, Feature
from typing import List, Iterator
import more_itertools
import pytest
 

def contains_at_least_n_features_n_bp_apart(feature_list: List[str], feature_count: int, distance_bp: int, count_multiple_copies: bool = False):
    """ 
    The operon must have at least feature_count given features at most distance_bp apart. 
    In other words, checks if there exists at least one sub-operon that contains some 
    combination of the features in the list and requires that those features be
    at most distance_bp apart.
    """
    serialized_list = "|".join(feature_list)
    custom_repr = f'contains-at-least-n-features-n-bp-apart:{serialized_list}-{feature_count}-{distance_bp}-{count_multiple_copies}'
    rule = Rule('contains_at_least_n_features_n_bp_apart',
                _contains_at_least_n_features_n_bp_apart,
                feature_list,
                feature_count,
                distance_bp,
                count_multiple_copies,
                custom_repr=custom_repr)
    return rule
    

def contains_these_features_n_bp_apart(feature_list: List[str], distance_bp: int, allow_gaps: bool = False):
    """
    All features in the list must appear grouped together in at least one location in the operon. 
    `distance_bp` sets the maximum allowed distance apart. By default, the features must be contiguous; 
    setting `allow_gaps=True` permits other annotated features to interupt the sequence if they also 
    adhere to the `distance_bp` requirement.
    """
    serialized_list= "|".join(feature_list)
    custom_repr = f'contains-these_features_n_bp_apart:{serialized_list}-{distance_bp}-{allow_gaps}'
    rule = Rule('contains-these_features_n_bp_apart',
                _contains_these_features_n_bp_apart,
                feature_list,
                distance_bp,
                allow_gaps,
                custom_repr=custom_repr)
    return rule
    

def _contains_at_least_n_features(operon: Operon, feature_names: List[str], feature_count: int, count_multiple_copies: bool) -> bool:
    """ Whether the operon has at least feature_count given features. """
    matches = [feature_name for feature_name in operon.feature_names if feature_name in feature_names]
    if len(matches) >= feature_count and count_multiple_copies:
        return True
    elif len(set(matches)) >= feature_count:
        return True
    else:
        return False
    

def _contains_at_least_n_features_n_bp_apart(operon: Operon, feature_list: List[str], feature_count: int, distance_bp: int, count_multiple_copies: bool) -> bool:
    """ 
    Whether the operon has at least feature_count given features at most distance_bp apart. 
    """
    # For this rule to be true, contains_at_least_n_features must be true also
    if not _contains_at_least_n_features(operon, feature_list, feature_count, count_multiple_copies):
        return False
    for sub_operon in _feature_clusters_with_same_orientation(operon, distance_bp):
        features_in_list = [feature for feature in sub_operon if feature in feature_list]
        if (len(features_in_list) >= feature_count and count_multiple_copies) or len(set(features_in_list)) >= feature_count:
            return True
    return False


def _contains_these_features_n_bp_apart(operon: Operon, feature_list: List[str], distance_bp: int, allow_gaps: bool) -> bool:
    """
    Whether all features in the list appear grouped together in at least one location in the operon. 
    """
    # For this rule to be true, contains_at_least_n_features must be true also
    if not _contains_at_least_n_features(operon, feature_list, feature_count=len(feature_list), count_multiple_copies=False):
        return False
    for sub_operon in _feature_clusters_with_same_orientation(operon, distance_bp):
        if set(sub_operon) == set(feature_list) or (set(feature_list).issubset(sub_operon) and allow_gaps):
            return True
    return False


def _same_orientation_two_features(f1: Feature, f2: Feature) -> bool:
    """ Checks if two features are in the same orientation. """
    # Somewhat arbitrarily, we say that two features cannot have the same 
    # orientation if one is a CRISPR array or some other directionless entity
    if (f1.strand is None) or (f2.strand is None):
        return False
    return (f1.strand == f2.strand)


def _feature_clusters_with_same_orientation(operon: Operon, distance_bp: int) -> Iterator[List[str]]:
    """
    Generates lists of names of features that appear in clusters in the
    operon. For a feature to be added to a cluster, it must be in the same 
    orientation as the cluster and at most distance_bp from its nearest
    upstream neighbor.
    """
    operon_features = [feature for feature in operon]
    operon_features.sort(key = lambda feature: feature.start)
    last_feature = operon_features.pop(0)
    cluster = [last_feature.name]
    for this_feature in operon_features:
        if _feature_distance(this_feature, last_feature) <= distance_bp and _same_orientation_two_features(this_feature, last_feature):
            cluster.append(this_feature.name)
        else:
            yield cluster
            cluster = [this_feature.name]
        last_feature = this_feature
    yield cluster


def contains_subset_group(feature_names: List[str], max_gap_distance_bp: int, require_same_orientation: bool, window_size: int):
    """ 
    Determines whether Features with the names in feature_names occur in any order without 
    interruption by other Features, with gaps between each Feature no larger than 
    max_gap_distance_bp. Groups do not need to contain all Features in feature_names, but the
    set of Features in the group must be <= the set of Features in feature_names. 
    """
    rule = Rule('contains_subset_group',
                _contains_subset_group,
                feature_names,
                max_gap_distance_bp,
                require_same_orientation,
                window_size)
    return rule


def _contains_subset_group(operon: Operon, feature_names: List[str], max_gap_distance_bp: int, require_same_orientation: bool, window_size: int) -> bool:
    """ 
    Determines whether Features with the names in feature_names occur in any order without 
    interruption by other Features, with gaps between each Feature no larger than 
    max_gap_distance_bp. Groups do not need to contain all Features in feature_names, but the
    set of Features in the group must be <= the set of Features in feature_names. 
    """
    assert window_size > 1
    if len(operon) < window_size:
        return False

    for operon_chunk in more_itertools.windowed(operon, window_size):
        operon_chunk_names = [feature.name for feature in operon_chunk]
        if set(operon_chunk_names) <= set(feature_names):
            max_gap_distance_in_group = 0
            for feature1, feature2 in zip(operon_chunk, operon_chunk[1:]):
                max_gap_distance_in_group = max(_feature_distance(feature1, feature2), max_gap_distance_in_group)
            if max_gap_distance_in_group > max_gap_distance_bp:
                continue
            if require_same_orientation:
                strands = set((feature.strand for feature in operon_chunk if feature.strand is not None))
                if len(strands) != 1:
                    continue
            return True
    return False


# TESTS

@pytest.mark.parametrize('gene1_start,gene1_end,gene2_start,gene2_end,strand1,strand2,expected', [
    (12, 400, 410, 600, 1, 1, True),
    (410, 600, 610, 400, 1, -1, False),
    (400, 500, 510, 600, 1, None, False),
    (400, 500, 510, 600, None, None, False)
    ])
def test_same_orientation_two_features(gene1_start, gene1_end, gene2_start, gene2_end, strand1, strand2, expected):
    f1 = Feature('f1', (gene1_start, gene1_end), '', strand1, '', 4e-19, 'a good gene', 'MCGYVER')
    f2 = Feature('f2', (gene2_start, gene2_end), '', strand2, '', 2e-5, 'a good gene', 'MGFRERAR')
    result = _same_orientation_two_features(f1, f2)
    assert result == expected


@pytest.mark.parametrize('gene_list,feature_count,distance_bp,count_multiple,expected', [
    (["cas1", "cas2", "cas3", "cas4"], 4, 100, False, True),
    (["cas1", "cas2", "cas3", "cas4"], 5, 100, False, False),
    (["cas1", "cas2", "cas3", "cas4"], 5, 100, True, True),
    (["cas1", "cas2", "cas3", "cas4"], 5, 10, True, False),
    (["cas1", "cas2", "cas3", "cas4", "cas5"], 5, 100, False, False)
    ])
def test_contains_at_least_n_features_n_bp_apart(gene_list, feature_count, distance_bp, count_multiple, expected):
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (450, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas3', (650, 680), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (700, 730), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas4', (740, 800), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas5', (900, 810), 'lcl|410|600|1|-1', -1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')
            ]
    operon = Operon('QCDRTU', "a file", 0, 3400, genes)
    rs = RuleSet().custom(rule=contains_at_least_n_features_n_bp_apart(gene_list, feature_count, distance_bp, count_multiple))
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize('gene_list,distance_bp,allow_gaps,expected', [
    (["tniQ", "cas1", "cas2", "cas4"], 100, True, True),
    (["tniQ", "cas1", "cas2", "cas3"], 100, False, False),
    (["tniQ", "cas4"], 100, False, False)
    ])
def test_contains_these_features_n_bp_apart(gene_list, distance_bp, allow_gaps, expected):
    genes = [
            Feature('tnsA', (20, 80), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('tniQ', (200, 340), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas1', (350, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (450, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas3', (650, 680), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (700, 730), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas4', (740, 800), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas5', (900, 810), 'lcl|410|600|1|-1', -1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')
            ]
    operon = Operon('QCDRTU', "a file", 0, 3400, genes)
    rs = RuleSet().custom(rule=contains_these_features_n_bp_apart(gene_list, distance_bp, allow_gaps))
    result = rs.evaluate(operon)
    assert result.is_passing is expected
