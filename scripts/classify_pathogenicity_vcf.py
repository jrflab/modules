def is_fs_splice_stop(record):
    ann_effect = get_ann_effect(record)
    return any([c in ef for ef in ann_effect for c in ["frameshift", "splice_donor", "splice_acceptor", "stop_gained"]])

def is_synonymous(record):
    ann_effect = get_ann_effect(record)
    return any(["synonymous_variant" in ef for ef in ann_effect])


def is_loh(record):
    return 'facetsLOH' in record.INFO


def is_provean_pathogenic(record):
    return 'provean_pred' in record.INFO and record.INFO['provean_pred'] == 'Deleterious'


def is_mt_pathogenic(record):
    return 'MT_pred' in record.INFO and 'disease' in record.INFO['MT_pred']

def is_dbnsfp_mt_passenger(record):
    return 'dbNSFP_MutationTaster_pred' not in record.INFO or 'P' in record.INFO['dbNSFP_MutationTaster_pred'] or \
        'N' in record.INFO['dbNSFP_MutationTaster_pred']


def is_chasm_pathogenic(record):
    cp = filter(lambda x: x.endswith('chasm_score'), record.INFO.keys())
    chasm_scores = [min(record.INFO[x]) for x in cp]
    return any([x <= 0.3 for x in chasm_scores])


def is_cancer_gene(record):
    return 'lawrence' in record.INFO or 'kandoth' in record.INFO or 'cancer_gene_census' in record.INFO


def get_ann_effect(record):
    return [x.split('|')[1] for x in record.INFO['ANN']]


def is_inframe(record):
    ann_effect = get_ann_effect(record)
    return any(["inframe" in ef for ef in ann_effect])


def is_missense(record):
    ann_effect = get_ann_effect(record)
    return any(["missense_variant" in ef for ef in ann_effect])


def is_hotspot(record):
    return "HOTSPOT" in record.INFO or 'hotspot' in record.INFO


def is_hap_insuf(record):
    return 'hap_insuf' in record.INFO


def is_fathmm_pathogenic(record):
    return "CANCER" in record.INFO['fathmm_pred'] if 'fathmm_pred' in record.INFO else False


def get_fs_splice_stop_pathogenicity(record):
    if (is_loh(record) or is_hap_insuf(record)) or is_cancer_gene(record):
        return "likely_pathogenic"
    else:
        return "passenger"


def get_missense_pathogenicity(record):
    if is_chasm_pathogenic(record):
        return "likely_pathogenic"
    elif is_dbnsfp_mt_passenger(record):
        return "passenger"
    elif is_fathmm_pathogenic(record) or is_hotspot(record):
        return "likely_pathogenic"
    else:
        return "passenger"


def get_inframe_pathogenicity(record):
    if (is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record)) and \
            is_mt_pathogenic(record) and is_provean_pathogenic(record):
        return 'likely_pathogenic'
    else:
        return 'passenger'


def classify_pathogenicity(record):
    """ classify pathogenicity for a vcf record
    """
    if is_missense(record):
        record.INFO['pathogenicity'] = get_missense_pathogenicity(record)
    elif is_fs_splice_stop(record):
        record.INFO["pathogenicity"] = get_fs_splice_stop_pathogenicity(record)
    elif is_inframe(record):
        record.INFO["pathogenicity"] = get_inframe_pathogenicity(record)
    elif is_synonymous(record):
        record.INFO["pathogenicity"] = 'passenger'
    else:
        record.INFO["pathogenicity"] = None


def requires_mt_provean(record):
    return is_inframe(record) and (is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record))
