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
parser.add_argument('--pass_only', action='store_true', help='output only somatic variants', default=False)
parser.add_argument('vcf_infile')

args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

vcf_reader.filters['germline5x'] = vcf.parser._Filter(id='germline5x',
                                                     desc='filter if normal depth > 20 '
                                                     'and normal VAF > 1/5 * tumor VAF '
                                                     'or normal variant depth greater than 1')
vcf_reader.filters['depth'] = vcf.parser._Filter(id='depth',
                                                 desc='filter if normal depth < 5 or tumor depth < 5')
vcf_reader.filters['noAF'] = vcf.parser._Filter(id='noAF', desc='missing AF')

vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

for record in vcf_reader:
    if record.FILTER is None:
        record.FILTER = []
    if record.genotype(args.tumor).data.AD is None or record.genotype(args.normal).data.AD is None or \
            sum(record.genotype(args.normal).data.AD) == 0 or sum(record.genotype(args.tumor).data.AD) == 0:
        record.FILTER.append('noAF')
    else:
        assert len(record.genotype(args.normal).data.AD) >= 2
        normal_ao = record.genotype(args.normal).data.AD[1:]
        normal_dp = sum(record.genotype(args.normal).data.AD)
        tumor_ao = record.genotype(args.tumor).data.AD[1:]
        tumor_dp = sum(record.genotype(args.tumor).data.AD)
        if tumor_dp < 5 or normal_dp < 5:
            record.FILTER.append('depth')
        if normal_dp != 0 and tumor_dp != 0:
            max_tumor_af = max([float(x) / tumor_dp for x in tumor_ao])
            max_normal_af = max([float(x) / normal_dp for x in normal_ao])
            if normal_dp > 20:
                if max_normal_af > float(max_tumor_af) / 5:
                    record.FILTER.append('germline5x')
            elif max(normal_ao) > 1:
                record.FILTER.append('germline5x')
        else:
            record.FILTER.append('noAF')
    if not args.pass_only or len(record.FILTER) == 0:
        vcf_writer.write_record(record)

vcf_writer.close()
