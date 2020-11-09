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


_no_scarlet_mutation_filter = rules.Filter('no-scarlet', _no_scarlet)
_arrays_inside_cds_filter = rules.Filter('array-is-inside-cds', _remove_arrays_inside_cds)

fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9) \
                      .must_be_within_n_bp_of_anything(30000) \
                      .custom(_arrays_inside_cds_filter) \
                      .custom(_no_scarlet_mutation_filter)
