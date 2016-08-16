#!/usr/bin/env python

""" split_snps_indels_vcf.py
split a vcf file into snps and indels/everything else
"""

import argparse
import vcf

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='split_snps_indels_vcf.py',
                                     description='split vcf file into snps and indels')
    parser.add_argument('vcf_infile')
    parser.add_argument('--snps', '-s', nargs='?', required=True, help='snp output vcf file')
    parser.add_argument('--indels', '-i', nargs='?', required=True, help='indel/everything else output vcf file')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    snp_vcf_writer = vcf.Writer(open(args.snps, 'w'), vcf_reader)
    indel_vcf_writer = vcf.Writer(open(args.indels, 'w'), vcf_reader)

    for record in vcf_reader:
        if record.is_snp:
            snp_vcf_writer.write_record(record)
        else:
            indel_vcf_writer.write_record(record)

    snp_vcf_writer.close()
    indel_vcf_writer.close()
