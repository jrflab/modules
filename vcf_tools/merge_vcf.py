#!/usr/bin/env python

import argparse
import vcf
import sys
from itertools import chain
from itertools import izip_longest

parser = argparse.ArgumentParser(prog='merge_vcf.py',
                                 description='merge vcf files')
parser.add_argument('vcf_files', nargs='+', help='vcf files to merge')
parser.add_argument('--out_file', nargs='?', help='merged vcf output file', default=sys.stdout)

args = parser.parse_args()

vcf_readers = [vcf.Reader(open(f, 'r')) for f in args.vcf_files]
vcf_reader = vcf_readers[0]
# merge header
if len(vcf_readers) > 1:
    for vcf_reader2 in vcf_readers[1:]:
        for inf in vcf_reader2.infos:
            if inf not in vcf_reader.infos:
                vcf_reader.infos[inf] = vcf_reader2.infos[inf]
        for filt in vcf_reader2.filters:
            if filt not in vcf_reader.infos:
                vcf_reader.filters[filt] = vcf_reader2.filters[filt]

vcf_writer = vcf.Writer(open(args.out_file, 'w'), vcf_reader)

sentinel = object()
for records in izip_longest(*vcf_readers, fillvalue=sentinel):
    if sentinel in records:
        raise ValueError('vcf files have different lengths')
    record = records[0]
    if len(records) > 1:
        #import pdb; pdb.set_trace()
        print [r.CHROM + str(r.POS) for r in records]
        assert all([r.CHROM == record.CHROM and r.POS == record.POS for r in records[1:]])
        ids = [r.ID for r in records]
        ids = [x.strip() for x in filter(None, ids)]
        ids = set(chain.from_iterable([i.split(';') for i in ids]))
        ids = filter(lambda x: x != '.', ids)
        record.ID = ';'.join(ids) if len(ids) > 0 else '.'
        for record2 in records[1:]:
            # merge filters in
            for f in record2.FILTER:
                if f not in record.FILTER:
                    record.add_filter(f)
            # merge info fields in
            for inf in record2.INFO.keys():
                if inf not in record.INFO:
                    record.add_info(inf, record2.INFO[inf])
    if len(record.FILTER) == 0:
        vcf_writer.write_record(record)

vcf_writer.close()
