from operon_analyzer import rules, genes


def _no_scarlet(operon: genes.Operon, ignored_reason_message: str):
    """ There is an ABC transporter gene that was accidentally included in our database because its truncation was
    created with Cas9, and Cas9 was in its gene name. """
    for feature in operon.get('cas9'):
        if 'ABC' in feature.description:
            feature.ignore(ignored_reason_message)


def _remove_arrays_inside_cds(operon: genes.Operon, ignored_reason_message: str):
    """
    Filters out CRISPR arrays that are contained within a CDS. These occur frequently
    due to the enthusiasm of pilercr.
    """
    genes = [feature for feature in operon.all_features if feature.name != 'CRISPR array']
    for array in operon:
        if array.name != 'CRISPR array':
            continue
        for gene in genes:
            if rules._feature_distance(gene, array) == 0:
                array.ignore(ignored_reason_message)
                break


def _pick_tns_protein_if_overlapped(operon: genes.Operon, ignored_message: str):
    """
    If a gene BLASTs as both a Tn7 protein and anything else, we prefer the Tn7 annotation,
    even if it has a lower e-value. This is because extremely divergent genes (of which
    there are many) will not have high e-values. This is still a heuristic, and we must
    later confirm that the identity is correct via phylogenetic analysis.
    """
    tn7_genes = ('tnsB', 'tnsC', 'tniQ', 'cas12k')
    for feature in operon.all_genes:
        for other_feature in operon.all_genes:
            if feature is other_feature:
                # don't compare feature to itself
                continue
            overlap = rules._calculate_overlap(feature, other_feature)
            if overlap is None:
                # these features do not overlap
                continue
            if overlap >= 0.9:
                if other_feature.name in tn7_genes and feature.name not in tn7_genes:
                    feature.ignore("Tn7 gene ignored")
                elif feature.bit_score > other_feature.bit_score:
                    other_feature.ignore("Tn7 gene ignored")


_no_scarlet_mutation_filter = rules.Filter('no-scarlet', _no_scarlet)
_arrays_inside_cds_filter = rules.Filter('array-is-inside-cds', _remove_arrays_inside_cds)
_pick_tns_filter = rules.Filter('pick-tn7', _pick_tns_protein_if_overlapped)


fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9) \
                      .must_be_within_n_bp_of_anything(30000) \
                      .custom(_arrays_inside_cds_filter) \
                      .custom(_no_scarlet_mutation_filter)

tn7fs = rules.FilterSet().must_be_within_n_bp_of_anything(30000) \
                      .custom(_arrays_inside_cds_filter) \
                      .custom(_no_scarlet_mutation_filter) \
                      .custom(_pick_tns_filter)
