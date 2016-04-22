#!/usr/bin/env python
# retain only novel variants

import argparse
import vcf
import re

parser = argparse.ArgumentParser(prog='common_filter_vcf.py',
                                 description='filter vcf file for novel variants (no rs-id/GMAF < 0.01 or cosmic id)')
parser.add_argument('vcf_infile')
parser.add_argument('vcf_outfile')

args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

vcf_reader.filters['Common'] = vcf.parser._Filter(id='common',
                                                  desc='no cosmic id, or has dbsnp ID and GMAF is > 0.01')

vcf_writer = vcf.Writer(open(args.vcf_outfile, 'w'), vcf_reader)

for record in vcf_reader:
    # ignore entries with cosmic IDs
    cosm_match = re.search(r'COSM', record.ID)
    if cosm_match is None:
        # filter entries with dbsnp IDs unless GMAF > 0.01
        rs_match = re.search(r'rs', record.ID)
        if rs_match is not None and ('GMAF' not in record.INFO or record.INFO['GMAF'] > 0.01):
            record.FILTER.append('Common')
    vcf_writer.write_record(record)

vcf_writer.close()
