#!/usr/bin/env python
""" add ups-coordinate to INFO of vcf file
"""

import argparse
import vcf
import pandas as pd
import sys
import copy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='merge_uvcf_vcf.py',
                                     description='add ups-coordinate to INFO of vcf file')
    parser.add_argument('uvcf_infile')
    parser.add_argument('vcf_infile')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))
    uvcf = pd.read_csv(args.uvcf_infile, comment='#', sep='\t')
    uvcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'UPS-COORDINATE', 'INFO']

    vcf_reader.infos['UPS_Coord'] = vcf.parser._Info(id='UPS_Coord', num='.', type='String',
                                                     desc="UPS-coordinate", source=None, version=None)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    ups_map = {}
    for i, row in uvcf.iterrows():
        x = '{}:{}_{}/{}'.format(row['CHROM'], row['POS'], row['REF'], row['ALT'])
        ups_map[x] = row['UPS-COORDINATE'].replace(" ", "")

    for record in vcf_reader:
        ups_coords = []
        for alt in record.ALT:
            alt = str(alt)
            x = '{}:{}_{}/{}'.format(record.CHROM, record.POS, record.REF, alt)
            if x in ups_map:
                ups_coords.append(ups_map[x])
            else:
                ups_coords.append('N/A[]')
        record.INFO['UPS_Coord'] = ups_coords
        vcf_writer.write_record(record)
    vcf_writer.close()
