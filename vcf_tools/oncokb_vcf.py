#!/usr/bin/env python

import argparse
import vcf
import sys
import pandas as pd
import re
try:
    from itertools import izip as zip
except ImportError:
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='oncokb_vcf.py',
                                     description='Add oncoKB annotation to vcf')
    parser.add_argument('--oncokb', help='oncoKB annotation file')
    parser.add_argument('vcf_infile')

    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    oncokb = pd.read_table(args.oncokb)

    vcf_reader.infos['oncoKB_level'] = vcf.parser._Info(id='oncoKB_level', num='.', type='String',
                                                        desc="OncoKB level(s)", source=None, version=None)
    vcf_reader.infos['oncoKB_cancer_type'] = vcf.parser._Info(id='oncoKB_cancer_type', num='.', type='String',
                                                              desc="OncoKB cancer type", source=None, version=None)

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    for record in vcf_reader:
        if 'SYMBOL' in record.INFO and 'HGVSp_Short' in record.INFO:
            assert len(record.INFO['SYMBOL']) == len(record.INFO['HGVSp_Short'])
            for symb, hgvsp in zip(record.INFO['SYMBOL'], record.INFO['HGVSp_Short']):
                hgvsp = re.sub(r'^p\.', r'', hgvsp)
                q = oncokb.query(
                    'Gene == "{}" and Alteration == "{}"'.format(symb, hgvsp))
                if len(q) > 0:
                    if 'oncoKB_level' not in record.INFO:
                        record.INFO['oncoKB_level'] = []
                        record.INFO['oncoKB_cancer_type'] = []
                    record.INFO['oncoKB_level'].extend(q['Level'])
                    record.INFO['oncoKB_cancer_type'].extend(map(lambda x: re.sub(' ', '_', x), q['Cancer Type']))
        vcf_writer.write_record(record)
    vcf_writer.close()
