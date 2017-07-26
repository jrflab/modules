#!/usr/bin/env python

import pandas as pd
import sys
import argparse
import math

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='compare_facets_cncf.py',
                                     description='compare two facets cncf files')
    parser.add_argument('cncf1')
    parser.add_argument('cncf2')
    args = parser.parse_args()

    df1 = pd.read_table(args.cncf1).fillna(0)
    df2 = pd.read_table(args.cncf2).fillna(0)
    mafr_diff = 0
    for i, row1 in df1.iterrows():
        for j, row2 in df2.iterrows():
            if (row1['loc.start'] >= row2['loc.start'] and
                    row1['loc.start'] <= row2['loc.end']) or \
                    (row1['loc.end'] >= row2['loc.start'] and
                     row1['loc.end'] <= row2['loc.start']):
                mafr_diff += math.fabs(row1['mafR.clust'] - row2['mafR.clust'])
                break
    if mafr_diff < 20:
        print(("success, CNCF files are similar: {} {}".format(args.cncf1, args.cncf2)))
        sys.exit(0)
    else:
        print(("failed, files have high mafR difference: {} {}".format(args.cncf1, args.cncf2)))
        sys.exit(1)
