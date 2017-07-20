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
    parser = argparse.ArgumentParser(prog='bed_annotate_vcf.py',
                                     description='annotate vcf file using a bed file')
    parser.add_argument('--info_tag', help='info tag for annotation')
    parser.add_argument('interval_file')
    parser.add_argument('vcf_infile')
    args = parser.parse_args()

    intervals = pd.read_table(args.interval_file, header=None, dtype={0: str, 1: np.int32, 2: np.int32})
    intervals = intervals.rename(columns={0: 'chr', 1: 'start', 2: 'end'})
    trees = {}
    for chrom, interval in intervals.groupby('chr'):
        chrom = re.sub(r'chr', '', chrom)
        trees[chrom] = intervaltree.IntervalTree.from_tuples(list(zip(interval.start, interval.end)))

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))
    vcf_reader.infos[args.info_tag] = vcf.parser._Info(id=args.info_tag, num='0', type='Flag',
                                                       desc='{}: overlap'.format(args.info_tag),
                                                       source=None,
                                                       version=None)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    for record in vcf_reader:
        chrom = re.sub(r'chr', '', record.CHROM)
        if record.FILTER is None:
            record.FILTER = []
        if chrom in trees:
            query = trees[chrom].search(record.POS)
            if len(query) == 0 & len(record.REF) > 1:
                query = trees[chrom].search(record.POS + len(record.REF))
            if len(query) != 0:
                record.INFO[args.info_tag] = True
        vcf_writer.write_record(record)

    vcf_writer.close()
