#!/usr/bin/env python

import argparse
import vcf
import sys

parser = argparse.ArgumentParser(prog='merge_vcf.py',
                                 description='merge vcf files')
parser.add_argument('vcf_files', nargs='+', help='vcf files to merge')
parser.add_argument('--out_file', nargs='?', help='merged vcf output file', default=sys.stdout,
                    type=argparse.FileType('w'))

args = parser.parse_args()

vcf_readers = [vcf.Reader(open(f, 'r')) for f in args.vcf_files]
vcf_reader = vcf_readers[0]
# merge header
if len(vcf_readers) > 1:
    for vcf_reader2 in vcf_readers[1:]:
        for form in vcf_reader2.formats:
            if form not in vcf_reader.formats:
                vcf_reader.formats[form] = vcf_reader2.formats[form]
        for inf in vcf_reader2.infos:
            if inf not in vcf_reader.infos:
                vcf_reader.infos[inf] = vcf_reader2.infos[inf]
        for filt in vcf_reader2.filters:
            if filt not in vcf_reader.infos:
                vcf_reader.filters[filt] = vcf_reader2.filters[filt]

vcf_writer = vcf.Writer(args.out_file, vcf_reader)

cpra = set()
for vcf_reader in vcf_readers:
    for record in vcf_reader:
        chr_pos_ref_alt = ":".join([record.CHROM, str(record.POS), record.REF, str(record.ALT)])
        if chr_pos_ref_alt not in cpra:
            vcf_writer.write_record(record)
            cpra.add(chr_pos_ref_alt)

vcf_writer.close()
