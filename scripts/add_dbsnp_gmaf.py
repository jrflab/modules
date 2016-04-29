#!/usr/bin/env python
# add INFO field for RARE variants, GMAF < 0.01

import argparse
import vcf
import re
import sys

parser = argparse.ArgumentParser(prog='add_dbsnp_gmaf.py',
                                 description='add GMAF field to dbsnp using CAF field')
parser.add_argument('vcf_infile')

args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

vcf_reader.infos['GMAF'] = vcf.parser._Info(id='GMAF', num=1, type='Float',
                                            desc="global minor allele frequency from 1000g",
                                            source=None, version=None)


vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

for record in vcf_reader:
    if 'CAF' in record.INFO:
        caf_str = [re.sub(r'^\.$', '0', re.sub(r'[\[\]]', '', x)) for x in record.INFO['CAF'] if x is not None]
        caf = sorted([float(x) for x in caf_str])
        gmaf = caf[len(caf) - 2]
        record.INFO['GMAF'] = gmaf
    vcf_writer.write_record(record)

vcf_writer.close()
