#!/usr/bin/env python

import argparse
import vcf
import sys
import copy
import collections
from itertools import chain
from itertools import izip_longest

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='merge_vcf.py',
                                     description='merge vcf files')
    parser.add_argument('vcf_files', nargs='+', help='vcf files to merge')
    parser.add_argument('--pass_only', action='store_true', default=False, help='write only passing variants')
    parser.add_argument('--ignore_filter', action='store_true', default=False, help='ignore filter column')
    parser.add_argument('--out_file', nargs='?', help='merged vcf output file', default=sys.stdout,
                        type=argparse.FileType('w'))

    args = parser.parse_args()

    vcf_readers = [vcf.Reader(open(f, 'r')) for f in args.vcf_files]
    vcf_reader = vcf_readers[0]
    # merge header
    if len(vcf_readers) > 1:
        for vcf_reader2 in vcf_readers:
            for inf in vcf_reader2.infos:
                if inf not in vcf_reader.infos:
                    vcf_reader.infos[inf] = vcf_reader2.infos[inf]
            for filt in vcf_reader2.filters:
                if filt not in vcf_reader.filters:
                    vcf_reader.filters[filt] = vcf_reader2.filters[filt]
            for fo in vcf_reader2.formats:
                if fo not in vcf_reader.formats:
                    vcf_reader.formats[fo] = vcf_reader2.formats[fo]
    vcf_reader.formats['AD'] = vcf.parser._Format(id='AD', num='.', type='Integer', desc='Allelic depth')

    vcf_writer = vcf.Writer(args.out_file, vcf_reader)
    pass_records = list()
    sentinel = object()
    for records in izip_longest(*vcf_readers, fillvalue=sentinel):
        if sentinel in records:
            raise ValueError('vcf files have different lengths')
        for rec in records:
            if rec.FILTER is None or args.ignore_filter or (hasattr(rec.FILTER, '__iter__') and
                                                            (len(rec.FILTER) == 0 or rec.FILTER[0] == '')):
                rec.FILTER = []
        record = records[0]
        new_record = copy.deepcopy(record)
        new_record.QUAL = None
        for sx in range(len(record.samples)):
            f_keys = record.FORMAT.split(':')
            f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
            if 'GT' not in f_keys or f_keys[0] != "GT":
                f_keys.insert(0, "GT")
                f_vals.insert(0, '0')
            elif f_vals[0] == '':
                f_vals[0] = '0'
            if 'AD' not in f_keys:
                if 'NR' in f_keys and 'NV' in f_keys:
                    f_keys.append("AD")
                    if type(record.samples[sx]['NR']) is list:
                        f_vals.append(int(record.samples[sx]['NR'][0]))
                        f_vals.append(map(lambda x: int(x), record.samples[sx]['NV']))
                    else:
                        f_vals.append([int(record.samples[sx]['NR']),
                                       int(record.samples[sx]['NV'])])
                if 'TIR' in f_keys and 'TAR' in f_keys:
                    f_keys.append("AD")
                    f_vals.append([int(record.samples[sx]['DP']) - int(record.samples[sx]['TIR'][0]),
                                   int(record.samples[sx]['TIR'][0])])
            new_record.samples[sx].data = collections.namedtuple('CallData', f_keys)
            handy_dict = dict(zip(f_keys, f_vals))
            new_vals = [handy_dict[x] for x in f_keys]
            # finally set CallData
            new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
        new_record.FORMAT = ":".join(new_record.samples[0].data._fields)

        if args.ignore_filter:
            new_record.FILTER = []
        if new_record.INFO is None or '.' in new_record.INFO:
            new_record.INFO = {}
        if len(records) > 1:
            assert all([r.CHROM == record.CHROM and r.POS == record.POS for r in records])

            ids = [r.ID for r in records]
            ids = [x.strip() for x in filter(None, ids)]
            ids = set(chain.from_iterable([i.split(';') for i in ids]))
            ids = filter(lambda x: x != '.', ids)
            new_record.ID = ';'.join(ids) if len(ids) > 0 else '.'
            for record2 in records:
                # merge filters in
                for f in record2.FILTER:
                    if f not in new_record.FILTER:
                        new_record.add_filter(f)
                # merge info fields in
                for inf in record2.INFO.keys():
                    if inf is not None and inf != '.' and inf not in new_record.INFO:
                        new_record.add_info(inf, record2.INFO[inf])
        if not args.pass_only or len(new_record.FILTER) == 0:
            vcf_writer.write_record(new_record)

    vcf_writer.close()
