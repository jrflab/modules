def is_fs_splice_stop(record):
    return any([c in vc for vc in record.INFO['Variant_Classification'] for c in ["Frame_Shift", "Splice_Site",
                                                                                  "Nonsense", "Translation_Start_Site",
                                                                                  "Nonstop"]])


def is_synonymous(record):
    return any([c in vc for vc in record.INFO['Variant_Classification'] for c in ["Silent", "Splice_Region"]])


def is_loh(record):
    return 'facetsLOH' in record.INFO


def is_provean_pathogenic(record):
    return 'provean_pred' in record.INFO and record.INFO['provean_pred'] == 'Deleterious'


def is_provean_missing(record):
    return 'provean_pred' in record.INFO and record.INFO['provean_pred'] == 'none'


def is_mt_pathogenic(record):
    return 'MT_pred' in record.INFO and 'disease' in record.INFO['MT_pred']


def is_mt_missing(record):
    return 'MT_pred' in record.INFO and record.INFO['MT_pred'] == "none"


def is_mt_passenger(record):
    return 'MutationTaster_pred' not in record.INFO or 'P' in record.INFO['MutationTaster_pred'] or \
        'N' in record.INFO['MutationTaster_pred'] or None in record.INFO['MutationTaster_pred']


def is_chasm_pathogenic(record):
    cp = filter(lambda x: x.endswith('chasm_score'), record.INFO.keys())
    chasm_scores = [min(record.INFO[x]) for x in cp]
    return any([x <= 0.3 for x in chasm_scores])


def is_cancer_gene(record):
    return 'lawrence' in record.INFO or 'kandoth' in record.INFO or 'cancer_gene_census' in record.INFO


def get_ann_effect(record):
    return [x.split('|')[1] for x in record.INFO['ANN']]


def is_inframe(record):
    return any(["In_Frame" in vc for vc in record.INFO['Variant_Classification']])


def is_missense(record):
    return any(["Missense" in vc for vc in record.INFO['Variant_Classification']])


def is_hotspot(record):
    return "HOTSPOT" in record.INFO or 'hotspot' in record.INFO


def is_hap_insuf(record):
    return 'hap_insuf' in record.INFO


def is_fathmm_pathogenic(record):
    return "D" in record.INFO['FATHMM_pred'] if 'FATHMM_pred' in record.INFO else False


def get_fs_splice_stop_pathogenicity(record):
    if (is_loh(record) or is_hap_insuf(record)) or is_cancer_gene(record):
        return "likely_pathogenic"
    else:
        return "passenger"


def get_missense_pathogenicity(record):
    if not is_chasm_pathogenic(record) and is_mt_passenger(record):
        return "passenger"
    elif is_chasm_pathogenic(record) or is_fathmm_pathogenic(record) or is_hotspot(record):
        return "likely_pathogenic"
    else:
        return "passenger"


def get_inframe_pathogenicity(record, no_remote):
    if no_remote:
        if is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record):
            return 'potentially_pathogenic'
        else:
            return 'passenger'
    elif (is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record)) and \
            (is_mt_pathogenic(record) and is_provean_pathogenic(record)) or \
            (is_mt_missing(record) and is_provean_pathogenic(record)) or \
            (is_mt_pathogenic(record) and is_provean_missing(record)):
        return 'likely_pathogenic'
    else:
        return 'passenger'


def classify_pathogenicity(record, no_remote=False):
    """ classify pathogenicity for a vcf record
    """
    if is_missense(record):
        record.INFO['pathogenicity'] = get_missense_pathogenicity(record)
    elif is_fs_splice_stop(record):
        record.INFO["pathogenicity"] = get_fs_splice_stop_pathogenicity(record)
    elif is_inframe(record):
        record.INFO["pathogenicity"] = get_inframe_pathogenicity(record, no_remote)
    elif is_synonymous(record):
        record.INFO["pathogenicity"] = 'passenger'
    else:
        record.INFO["pathogenicity"] = None


def requires_mt_provean(record):
    return is_inframe(record) and (is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record))
