#!/usr/bin/env python

""" classify pathogenicity of vcf records, querying provean as necessary
"""

import argparse
import vcf
import sys


def get_ann_effect(record):
    return [x.split('|')[1] for x in record.INFO['ANN']]


def is_missense(record):
    ann_effect = get_ann_effect(record)
    return any(["missense_variant" in ef for ef in ann_effect])


def is_mt_passenger(record):
    return 'P' in record.INFO['dbNSFP_MutationTaster_pred'] or \
        'N' in record.INFO['dbNSFP_MutationTaster_pred']


def is_chasm_pathogenic(record):
    cp = filter(lambda x: x.endswith('chasm_score'), record.INFO.keys())
    chasm_scores = [min(record.INFO[x]) for x in cp]
    return any([x <= 0.3 for x in chasm_scores])


def is_hotspot(record):
    return "HOTSPOT" in record.INFO or 'hotspot' in record.INFO


def is_fathmm_pathogenic(record):
    return "CANCER" in record.INFO['fathmm_pred'] if 'fathmm_pred' in record.INFO else False


def is_hap_insuf(record):
    return 'hap_insuf' in record.INFO


def get_missense_pathogenicity(record):
    if is_chasm_pathogenic(record):
        return "likely_pathogenic"
    elif is_mt_passenger(record):
        return "passenger"
    elif is_fathmm_pathogenic(record) or is_hotspot(record):
        return "likely_pathogenic"
    else:
        return "passenger"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='classify_pathogenicity_vcf.py',
                                     description='Add pathogenicity to vcf file')
    parser.add_argument('vcf_infile')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    assert "ANN" in vcf_reader.infos
    assert "HOTSPOT" in vcf_reader.infos or "hotspot" in vcf_reader.infos
    assert "fathmm_pred" in vcf_reader.infos
    assert "dbNSFP_MutationTaster_pred" in vcf_reader.infos

    # add necessary info headers
    vcf_reader.infos['pathogenicity'] = vcf.parser._Info(id='pathogenicity', num=-1, type='String',
                                                         desc="Classification of pathogenicity",
                                                         source=None, version=None)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    records = list()
    for record in vcf_reader:
        if is_missense(record):
            record.INFO["pathogenicity"] = get_missense_pathogenicity(record)
        vcf_writer.write_record(record)
    vcf_writer.close()
