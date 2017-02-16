#!/usr/bin/env python
# retain only PASS variants

import argparse
import vcf
import re
import sys

parser = argparse.ArgumentParser(prog='pass_filter_vcf.py',
                                 description='filter vcf file for PASS variants')
parser.add_argument('vcf_infile')

args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

for record in vcf_reader:
    if len(record.FILTER) == 0:
        vcf_writer.write_record(record)

vcf_writer.close()
