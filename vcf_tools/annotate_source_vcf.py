#!/usr/bin/env python
""" add variant source (name of caller) to a vcf
"""

import argparse
import vcf
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='annotate_source_vcf.py',
                                     description='source annotation to add to each variant INFO')
    parser.add_argument('--source', required=True, help='name of source')
    parser.add_argument('--vcf_infile', required=False, type=argparse.FileType('r'), default=sys.stdin)

    args = parser.parse_args()

    vcf_reader = vcf.Reader(args.vcf_infile)

    vcf_reader.infos['variantCaller'] = vcf.parser._Info(id='variantCaller', num='.', type='String',
                                                         desc="variant caller(s) used to find the variant",
                                                         source=None, version=None)

    if args.source == 'lancet':
        vcf_reader.infos['LEN'] = vcf.parser._Info(id='LEN', num='1', type='Integer',
                                                   desc="length of insertion/deletion",
                                                   source=None, version=None)
        vcf_reader.infos['TYPE'] = vcf.parser._Info(id='TYPE', num='1', type='String',
                                                    desc="insertion or deletion",
                                                    source=None, version=None)

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    for record in vcf_reader:
        record.INFO['variantCaller'] = [args.source]
        vcf_writer.write_record(record)
    vcf_writer.close()
