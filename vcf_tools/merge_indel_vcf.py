#!/usr/bin/env python
""" merge multiple vcf files that have uvcf annotation (as UPS_Coord)
"""

import argparse
import vcf
import sys
from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='merge_indel_vcf.py',
                                     description='merge multiple indel uvcf files')
    parser.add_argument('vcf_file', nargs='+', help='indel vcf files to marge')
    parser.add_argument('--ignore_filter', action='store_true', default=False, help='ignore filter column')
    parser.add_argument('--min_keep', type=int, default=2, help='minimum presence to keep variant')
    args = parser.parse_args()

    vcf_readers = [vcf.Reader(open(f, 'r')) for f in args.vcf_file]
    vr0 = vcf_readers[0]
    #assert all(['UPS_Coord' in vr.infos for vr in vcf_readers])

    vr0.filters['minCaller'] = vcf.parser._Filter(id='minCaller',
                                                  desc='variant is not present in at least '
                                                  '{} callers'.format(args.min_keep))
    vr0.infos['variantCaller'] = vcf.parser._Info(id='variantCaller', num='.', type='String',
                                                  desc="variant caller(s) used to find the variant",
                                                  source=None, version=None)
    if len(vcf_readers) > 1:
        for vcf_reader2 in vcf_readers:
            for inf in vcf_reader2.infos:
                if inf not in vr0.infos:
                    vr0.infos[inf] = vcf_reader2.infos[inf]
            for filt in vcf_reader2.filters:
                if filt not in vr0.infos:
                    vr0.filters[filt] = vcf_reader2.filters[filt]
            for fo in vcf_reader2.formats:
                if fo not in vr0.formats:
                    vr0.formats[fo] = vcf_reader2.formats[fo]

    vr0.formats['AD'] = vcf.parser._Format(id='AD', num='.', type='Integer', desc='Allelic depth')

    ups_map = defaultdict(list)
    i = 0
    for vcf_reader in vcf_readers:
        ups = set()
        for rec in vcf_reader:
            if rec.FILTER is None or args.ignore_filter or (hasattr(rec.FILTER, '__iter__') and
                                                            (len(rec.FILTER) == 0 or rec.FILTER[0] == '')):
                rec.FILTER = []
            for ups_coord in rec.INFO['UPS_Coord']:
                if ups_coord != 'N/A[]' and ups_coord not in ups:
                    ups_map[ups_coord].append(rec)
                    ups.add(ups_coord)

    vcf_writer = vcf.Writer(sys.stdout, vr0)
    for recs in ups_map.values():
        pass_recs = list(filter(lambda r: len(r.FILTER) == 0, recs))
        if len(pass_recs) >= args.min_keep:
            if 'variantCaller' not in recs[0].INFO:
                recs[0].INFO['variantCaller'] = []
            if len(pass_recs) > 1:
                for pr in pass_recs[1:]:
                    if 'variantCaller' in pr.INFO:
                        for vc in pr.INFO['variantCaller']:
                            recs[0].INFO['variantCaller'].append(vc)
            vcf_writer.write_record(recs[0])
        else:
            rec = recs[0]
            rec.FILTER.append('minCaller')
            vcf_writer.write_record(rec)
    vcf_writer.close()
