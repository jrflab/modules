#!/usr/bin/env python
""" annotate off-target and filter off-target with low depth
"""

import argparse
import vcf
import sys
import pandas as pd
import numpy as np
import intervaltree
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='interval_depth_filter_vcf.py',
                                     description='filter off-target with low depth')
    parser.add_argument('--depth_threshold', default=50, help='min depth for off-target calls')
    parser.add_argument('interval_file')
    parser.add_argument('vcf_infile')
    args = parser.parse_args()

    intervals = pd.read_table(args.interval_file, header=None, dtype={0: str, 1: np.int32, 2: np.int32})
    intervals = intervals.rename(columns={0: 'chr', 1: 'start', 2: 'end'})
    trees = {}
    for chrom, interval in intervals.groupby('chr'):
        chrom = re.sub(r'chr', '', chrom)
        trees[chrom] = intervaltree.IntervalTree.from_tuples(list(zip(interval.start - 50, interval.end + 50)))

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))
    vcf_reader.filters['lowDepthOffTarget'] = vcf.parser._Filter(
        id='lowDepthOffTarget',
        desc='low depth (< {}) and no overlap with intervals'.format(args.depth_threshold))
    vcf_reader.infos['offTarget'] = vcf.parser._Info(id='offTarget', num='0', type='Flag',
                                                     desc='no overlap with target intervals',
                                                     source=None,
                                                     version=None)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    for record in vcf_reader:
        chrom = re.sub(r'chr', '', record.CHROM)
        if getattr(record.samples[0].data, 'DP', None) is not None:
            depths = [x['DP'] for x in record.samples]
        elif getattr(record.samples[0].data, 'AD', None) is not None:
            depths = [sum(x['AD']) for x in record.samples]
        elif getattr(record.samples[0].data, 'NR', None) is not None:
            depths = [x['NR'] for x in record.samples]
        else:
            raise ValueError('no depth or allelic depth available')

        if record.FILTER is None:
            record.FILTER = []
        if chrom not in trees:
            record.INFO['offTarget'] = True
            if any([x < args.depth_threshold for x in depths]):
                record.FILTER.append('lowDepthOffTarget')
        else:
            query = trees[chrom].search(record.POS)
            if len(query) == 0:
                record.INFO['offTarget'] = True
                if any([x < args.depth_threshold for x in depths]):
                    record.FILTER.append('lowDepthOffTarget')
        vcf_writer.write_record(record)

    vcf_writer.close()
