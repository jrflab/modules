#!/usr/bin/env python
# filter SNPs that are rare and have a global minor allele frequency less than 0.01

import argparse
import vcf
import re
import sys

parser = argparse.ArgumentParser(prog='filter_dbsnp_gmaf.py',
                                 description='filter dbsnp gmaf (CAF info)')
parser.add_argument('vcf_infile')

args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

for record in vcf_reader:
    if 'CAF' not in record.INFO:
        vcf_writer.write_record(record)
    else:
        caf1 = float(re.sub(r'[\[\]]', '', record.INFO['CAF'][0]))
        if caf1 < 0.99:
            vcf_writer.write_record(record)
vcf_writer.close()
