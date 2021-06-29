import sys
from typing import List

from operon_analyzer import load, genes

from scan import scan_for_repeats, DNASlice, RepeatPair


SPACER_SEARCH_SIZE = 200
TARGET_SEARCH_SIZE = 500
MIN_TARGET_SIZE = 10
MAX_TARGET_SIZE = 24
UPSTREAM = "upstream"
DOWNSTREAM = "downstream"


def create_search_regions(operon: genes.Operon):
    for feature, seed_location, direction in define_search_seed(operon):
        bounds_result = find_array_bounds(operon, feature, seed_location, direction)
        if not bounds_result:
            continue
        start, end = bounds_result
        array_start, array_end = min(start, end), max(start, end)
        sequence = load.load_sequence(operon)

        result = create_sts_search_regions(array_start, array_end, seed_location, direction, sequence)
        if result is None:
            continue
        source, target = result
        operon._features.append(source)
        operon._features.append(target)


def find_self_targeting_spacers(operon: genes.Operon):
    sources = []
    targets = []
    for feature in operon:
        if feature.name == 'STS source':
            sources.append(feature)
        elif feature.name == 'STS target':
            targets.append(feature)

    sequence = load.load_sequence(operon)
    for source in sources:
        direction = source.description  # we overloaded the description with the orientation of the nearest Cas12k
        repeats = []
        for alignment_target in targets:
            for inverted in (False, True):
                source_slice = DNASlice(sequence, source.start, source.end)
                target_slice = DNASlice(sequence, alignment_target.start, alignment_target.end)
                repeats.extend(scan_for_repeats(source_slice, target_slice, inverted, MIN_TARGET_SIZE, MAX_TARGET_SIZE, 0.9))
        if repeats:
            repeats.sort(key=lambda r: -r.ar.score)
            best = repeats[0]
            spacer = best.slice1
            target = best.slice2
            repeat_start, repeat_end = (spacer.start - 20, spacer.start) if direction == DOWNSTREAM else (spacer.end, spacer.end + 20)
            repeat = sequence[repeat_start: repeat_end] if direction == DOWNSTREAM else str(sequence[repeat_start: repeat_end].reverse_complement())
            spacer_sequence = str(spacer.sequence) if direction == DOWNSTREAM else str(spacer.sequence.reverse_complement())
            self_targeting_repeat = genes.Feature('STR', (repeat_start, repeat_end), '', 1 if direction == DOWNSTREAM else -1, '', None, '', repeat)
            self_targeting_spacer = genes.Feature('STS', (spacer.start, spacer.end), '', 1 if direction == DOWNSTREAM else -1, '', None, f'for {target.start - spacer.start}', spacer_sequence)
            target_is_inverted = (direction == DOWNSTREAM and best.ar.inverted) or (direction == UPSTREAM and not best.ar.inverted)

            if not target_is_inverted:
                pam_start, pam_end = target.start - 6, target.start
                pam_sequence = sequence[pam_start: pam_end]
            else:
                pam_start, pam_end = target.end, target.end + 6
                pam_sequence = str(sequence[pam_start: pam_end].reverse_complement())
            self_targeting_pam = genes.Feature('PAM', (pam_start, pam_end), '', 1 if not target_is_inverted else -1, '', None, '', pam_sequence)

            self_target_sequence = str(target.sequence) if not target_is_inverted else str(target.sequence.reverse_complement())
            self_target = genes.Feature('ST', (target.start, target.end), '', 1 if not target_is_inverted else -1, '', None, f'from {spacer.start - target.start}', self_target_sequence)

            operon._features.append(self_targeting_repeat)
            operon._features.append(self_targeting_spacer)
            operon._features.append(self_targeting_pam)
            operon._features.append(self_target)


def create_sts_search_regions(array_start, array_end, cas12_boundary, direction, sequence) -> List[RepeatPair]:
    assert array_start < array_end
    if direction == UPSTREAM:
        source = DNASlice(sequence,
                          array_start - SPACER_SEARCH_SIZE,
                          cas12_boundary)
        target = DNASlice(sequence,
                          max(0, array_start - SPACER_SEARCH_SIZE - TARGET_SEARCH_SIZE),
                          array_start - SPACER_SEARCH_SIZE)

    elif direction == DOWNSTREAM:
        source = DNASlice(sequence,
                          cas12_boundary,
                          array_end + SPACER_SEARCH_SIZE)
        target = DNASlice(sequence,
                          array_end + SPACER_SEARCH_SIZE,
                          min(len(sequence), array_end + SPACER_SEARCH_SIZE + TARGET_SEARCH_SIZE))
    if len(source.sequence) < MAX_TARGET_SIZE or len(target.sequence) < MAX_TARGET_SIZE:
        return None
    # if the CRISPR array is too far from Cas12k, we assume it's not used
    if len(source.sequence) > 5000:
        return None

    source_feature = genes.Feature('STS source', (source.start, source.end), '', None, '', None, direction, source.sequence)
    target_feature = genes.Feature('STS target', (target.start, target.end), '', None, '', None, '', target.sequence)
    return source_feature, target_feature


def define_search_seed(operon):
    """ Finds Cas12k(s) and determines where to start searching. We start immediately downstream of the gene. """
    for feature in operon:
        if feature.name == 'cas12k':
            yield (feature, feature.end, DOWNSTREAM) if feature.strand == 1 else (feature, feature.start, UPSTREAM)


def find_array_bounds(operon: genes.Operon, cas12k: genes.Feature, cas12k_end: int, direction: str):
    assert direction in (DOWNSTREAM, UPSTREAM)
    arrays = []
    for feature in operon:
        is_array = feature.name in ('CRISPR array', 'Repeat Spacer')
        array_on_correct_side = ((feature.start > cas12k_end and direction == DOWNSTREAM) \
                                or (feature.start < cas12k_end and direction == UPSTREAM))
        if is_array and array_on_correct_side:
            arrays.append(abs(cas12k_end - feature.start))
            arrays.append(abs(cas12k_end - feature.end))
    if not arrays:
        return None
    arrays.sort()
    start = arrays[0]
    end = arrays[1]
    if len(arrays) == 2:
        return (cas12k_end - start, cas12k_end - end) if direction == UPSTREAM else (cas12k_end + start, cas12k_end + end)
    for a, b in zip(arrays[1:], arrays[2:]):
        if b - a < 250:
            end = b
        else:
            break
    return (cas12k_end - start, cas12k_end - end) if direction == UPSTREAM else (cas12k_end + start, cas12k_end + end)


if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        operons = list(load.load_operons(f))
        for operon in operons:
            create_search_regions(operon)
            find_self_targeting_spacers(operon)
            print(operon.as_str())
