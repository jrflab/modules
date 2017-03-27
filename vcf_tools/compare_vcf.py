#!/usr/bin/env python
"""
compare two vcf files after melting
"""

import pandas as pd
import numpy as np
import subprocess
import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-q', help='exit code 0 on no diff', default=False, action='store_true')
    parser.add_argument('--ignore_info', help='ignore info columns', default=False, action='store_true')
    parser.add_argument('--ignore_cols', help='ignore specified columns', nargs='*', default=None)
    parser.add_argument('vcf_file1')
    parser.add_argument('vcf_file2')
    args = parser.parse_args()
    if args.vcf_file1.endswith('gz'):
        cmd1 = 'zcat {} | vcf_melt'
    else:
        cmd1 = 'cat {} | vcf_melt'
    if args.vcf_file2.endswith('gz'):
        cmd2 = 'zcat {} | vcf_melt'
    else:
        cmd2 = 'cat {} | vcf_melt'
    melt1 = subprocess.Popen(cmd1.format(args.vcf_file1), stdout=subprocess.PIPE, shell=True)
    melt2 = subprocess.Popen(cmd2.format(args.vcf_file2), stdout=subprocess.PIPE, shell=True)

    df1 = pd.read_table(melt1.stdout, sep='\t', low_memory=False, na_values=['None', '.'])
    df2 = pd.read_table(melt2.stdout, sep='\t', low_memory=False, na_values=['None', '.'])

    if not df1.equals(df2):
        ne_stacked = (df1.fillna('.') != df2.fillna('.')).stack()
        changed = ne_stacked[ne_stacked]
        changed.index.names = ['id', 'col']
        diff_locs = np.where(df1.fillna('.') != df2.fillna('.'))
        changed_from = df1.values[diff_locs]
        changed_to = df2.values[diff_locs]
        df = pd.DataFrame({'from': changed_from, 'to': changed_to}, index=changed.index)
        if args.ignore_info:
            df = df.reset_index().ix[~df.reset_index()['col'].str.startswith('info.'), :]
        if args.ignore_cols is not None:
            df = df.reset_index().ix[~df.reset_index()['col'].isin(args.ignore_cols), :]
        if args.q and len(df) > 0:
            exit(1)
        if len(df) > 0:
            df.to_csv(sys.stdout, sep='\t')
    exit(0)
