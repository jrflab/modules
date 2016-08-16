#!/usr/bin/env python
""" retain only somatic variants
"""

import argparse
import vcf
import sys

parser = argparse.ArgumentParser(prog='somatic_ad_filter_vcf.py',
                                 description='filter vcf file for somatic variants')
parser.add_argument('--tumor', '-t', nargs='?', help='tumor sample name', default='TUMOR')
parser.add_argument('--normal', '-n', nargs='?', help='normal sample name', default='NORMAL')
parser.add_argument('vcf_infile')

args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

vcf_reader.filters['somaticAlleleDepth'] = vcf.parser._Filter(id='somaticAlleleDepth',
                                                              desc='filter if normal depth > 20 '
                                                              'and normal VAF > 1/5 * tumor VAF '
                                                              'or normal variant depth greater than 1')
vcf_reader.filters['noAF'] = vcf.parser._Filter(id='noAF', desc='missing AF')

vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

for record in vcf_reader:
    assert record.genotype(args.tumor).data.AF is not None
    assert record.genotype(args.normal).data.AF is not None
    assert len(record.genotype(args.normal).data.AD) >= 2
    tumor_af = record.genotype(args.tumor).data.AF
    normal_af = record.genotype(args.normal).data.AF
    normal_ao = record.genotype(args.normal).data.AD[1:]
    normal_dp = sum(record.genotype(args.normal).data.AD)
    if tumor_af is not None and normal_af is not None:
        if isinstance(tumor_af, list):
            tumor_af = max(tumor_af)
        if isinstance(normal_af, list):
            normal_af = max(normal_af)
        normal_ao = max(normal_ao)

        if normal_dp > 20:
            if normal_af > tumor_af / 5:
                record.FILTER.append('somaticAlleleDepth')
        elif normal_ao > 1:
            record.FILTER.append('somaticAlleleDepth')
    else:
        record.FILTER.append('noAF')
    vcf_writer.write_record(record)

vcf_writer.close()
