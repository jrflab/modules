#!/usr/bin/env python
""" retain variants with varian reads in the tumor
"""

import argparse
import vcf
import sys

parser = argparse.ArgumentParser(prog='somatic_ad_filter_vcf.py',
                                 description='filter vcf file for somatic variants')
parser.add_argument('--tumor', '-t', nargs='?', help='tumor sample name', default='TUMOR')
parser.add_argument('--normal', '-n', nargs='?', help='normal sample name', default='NORMAL')
parser.add_argument('--pass_only', action='store_true', help='output only somatic variants', default=False)
parser.add_argument('--vcf_infile', default=sys.stdin, type=argparse.FileType('r'))

args = parser.parse_args()

vcf_reader = vcf.Reader(args.vcf_infile)

vcf_reader.filters['noTumorVariantReads'] = vcf.parser._Filter(id='noTumorVariantReads',
                                                               desc='no tumor variant reads')
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
        if normal_dp != 0 and tumor_dp != 0:
            max_tumor_af = max([float(x) / tumor_dp for x in tumor_ao])
            if max_tumor_af == 0:
                record.FILTER.append('noTumorVariantReads')
        else:
            record.FILTER.append('noAF')
    if not args.pass_only or len(record.FILTER) == 0:
        vcf_writer.write_record(record)

vcf_writer.close()
