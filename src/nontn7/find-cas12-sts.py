# do a pairwise search between the source (defined above) and 25kb upstream of cas12a and 25kb downstream of the array
# create Features of top hits and source


import sys

from operon_analyzer import load, genes
from scan import DNASlice, scan_for_repeats


for n, operon in enumerate(load.load_operons(sys.stdin)):
    sequence = load.load_sequence(operon)
    cas12_coords = None
    crispr_coords = []
    for feature in operon:
        if 'cpf1' in feature.name.lower() or 'cpf1' in feature.description.lower() and len(feature) > 3000:
            if cas12_coords is not None and cas12_coords != (feature.start, feature.end):
                continue
            cas12_coords = feature.start, feature.end
        elif feature.name == 'Repeat Spacer' or feature.name == 'CRISPR array':
            crispr_coords.append((feature.start, feature.end))
    if cas12_coords is None:
        continue
    best_array = 0
    best_distance = 10000000
    cas12_start, cas12_end = cas12_coords
    for n, (array_start, array_end) in enumerate(crispr_coords):
        distance = min(abs(cas12_start - array_start),
                       abs(cas12_start - array_end),
                       abs(cas12_end - array_start),
                       abs(cas12_end - array_end))
        if distance < best_distance:
            best_array = n
            best_distance = distance
    if not crispr_coords:
        continue
    array_start, array_end = crispr_coords[best_array]

    best_coords = sorted([cas12_start, cas12_end, array_start, array_end])

    # pick the region where we might find atypical spacers
    # buffer by 50 bp on either side to help diminish accidental capture of CRISPR repeats
    spacer_region_start, spacer_region_end = best_coords[1] + 50, best_coords[2] - 50
    if spacer_region_end <= spacer_region_start:
        continue
    # spacer_region_sequence = sequence[spacer_region_start: spacer_region_end]

    upstart = max(0, min(best_coords) - 25000)
    upend = min(best_coords)
    # target_up = sequence[upstart: upend]

    downstart = max(best_coords)
    downend = min(len(sequence), max(best_coords) + 25000)
    # target_down = sequence[downstart: downend]

    spacer_source = DNASlice(sequence, spacer_region_start, spacer_region_end)
    target_upstream = DNASlice(sequence, upstart, upend)
    target_downstream = DNASlice(sequence, downstart, downend)

    min_length = 16
    max_length = 28
    min_homology = 0.85

    up_results = []
    down_results = []
    for inverted in (True, False):
        up_results.extend(scan_for_repeats(spacer_source, target_upstream, inverted, min_length, max_length, min_homology))
        down_results.extend(scan_for_repeats(spacer_source, target_downstream, inverted, min_length, max_length, min_homology))

    if not up_results and not down_results:
        continue
    # pick the best alignment from each side
    if up_results:
        up_result = sorted(up_results, key=lambda x: -x.ar.score)[0]

        up_spacer_feature = genes.Feature(f'Self-targeting spacer for {up_result.slice2.start}',
                                          (up_result.slice1.start, up_result.slice1.end),
                                          '',
                                          None,
                                          '',
                                          None,
                                          f'{up_result.ar.aligned_needle_sequence}\n{up_result.ar.comparison_string}\n{up_result.ar.aligned_haystack_sequence}',
                                          up_result.ar.aligned_needle_sequence)

        up_target_feature = genes.Feature(f'Target from {up_result.slice1.start}',
                                          (up_result.slice2.start, up_result.slice2.end),
                                          '',
                                          None,
                                          '',
                                          None,
                                          f'{up_result.ar.aligned_needle_sequence}\n{up_result.ar.comparison_string}\n{up_result.ar.aligned_haystack_sequence}',
                                          up_result.ar.aligned_haystack_sequence)
        operon._features.append(up_spacer_feature)
        operon._features.append(up_target_feature)

    if down_results:
        down_result = sorted(down_results, key=lambda x: -x.ar.score)[0]
        down_spacer_feature = genes.Feature(f'Self-targeting spacer for {down_result.slice2.start}',
                                            (down_result.slice1.start, down_result.slice1.end),
                                            '',
                                            None,
                                            '',
                                            None,
                                            f'{down_result.ar.aligned_needle_sequence}\n{down_result.ar.comparison_string}\n{down_result.ar.aligned_haystack_sequence}',
                                            down_result.ar.aligned_needle_sequence)

        down_target_feature = genes.Feature(f'Target from {down_result.slice1.start}',
                                            (down_result.slice2.start, down_result.slice2.end),
                                            '',
                                            None,
                                            '',
                                            None,
                                            f'{down_result.ar.aligned_needle_sequence}\n{down_result.ar.comparison_string}\n{down_result.ar.aligned_haystack_sequence}',
                                            down_result.ar.aligned_haystack_sequence)
        operon._features.append(down_spacer_feature)
        operon._features.append(down_target_feature)
    print(operon.as_str())
