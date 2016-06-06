#!/usr/bin/env python
""" recenter base counts of normal in merged base count file and clip depths
maintain the proportion of ATCGs
"""

import argparse
import pandas as pd
import sys
import numpy as np


bases = ['A', 'C', 'G', 'T']


def scale(y, c=True, sc=True):
    x = y.copy()

    if c:
        x -= x.mean()
    if sc and c:
        x /= x.std()
    elif sc:
        x /= np.sqrt(x.pow(2).sum().div(x.count() - 1))
    return x


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='clip_base_count.py',
                                     description='threshold the tumor/normal depth of a merged base count file while'
                                     'maintaining the proportion of ATCGs')
    parser.add_argument('--normal_center', default=100, type=int, help='recenter normal depth')
    # parser.add_argument('--tumor_center', default=200, type=int, help='recenter tumor depth')
    parser.add_argument('--threshold', default=800, type=int, help='depth threshold')
    args = parser.parse_args()

    df = pd.read_csv(sys.stdin, sep='\t')
    df = df.fillna(0)

    # recenter normal depth at args.normal_center
    normal_dp_z = scale(df['NOR.DP'])
    normal_dp = args.normal_center + normal_dp_z * df['NOR.DP'].std()
    normal_dp = np.clip(normal_dp, 0, args.threshold).astype(int)
    tn = 'NOR'
    for base in bases:
        df[tn + '.' + base + 'p'] = (df[tn + '.' + base + 'p'] /
                                     df[tn + '.DP'].astype(float) * normal_dp).fillna(0).astype(int)
        df[tn + '.' + base + 'n'] = (df[tn + '.' + base + 'n'] /
                                     df[tn + '.DP'].astype(float) * normal_dp).fillna(0).astype(int)
    df[tn + '.DP'] = normal_dp

    # recenter tumor depth at args.tumor_center
    # tumor_dp_z = scale(df['TUM.DP'])
    # tumor_dp = args.tumor_center + tumor_dp_z * df['TUM.DP'].std()
    # tumor_dp = np.clip(tumor_dp, 0, args.threshold).astype(int)
    # tn = 'TUM'
    # for base in bases:
        # df[tn + '.' + base + 'p'] = (df[tn + '.' + base + 'p'] /
                                     # df[tn + '.DP'].astype(float) * tumor_dp).fillna(0).astype(int)
        # df[tn + '.' + base + 'n'] = (df[tn + '.' + base + 'n'] /
                                     # df[tn + '.DP'].astype(float) * tumor_dp).fillna(0).astype(int)
    # df[tn + '.DP'] = tumor_dp

    df.to_csv(sys.stdout, sep='\t', index=False)
