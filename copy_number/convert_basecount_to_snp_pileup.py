#!/usr/bin/env python

import pandas as pd
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='convert_basecount_to_snp_pileup.py',
                                     description='convert old basecount facets files to snp-pileup')
    parser.add_argument('basecount')
    args = parser.parse_args()

    bc = pd.read_table(args.basecount, dtype={'Chrom': str})
    sp = pd.DataFrame()
    sp['Chromosome'] = bc.Chrom
    sp['Position'] = bc.Pos
    sp['Ref'] = bc.Ref
    sp['Alt'] = bc.Alt
    for i, row in bc.iterrows():
        sp.ix[i, 'File1R'] = int(bc.ix[i, "NOR.{}p".format(sp.ix[i, 'Ref'])] +
                                 bc.ix[i, "NOR.{}n".format(sp.ix[i, 'Ref'])])
        sp.ix[i, 'File1A'] = int(bc.ix[i, "NOR.{}p".format(sp.ix[i, 'Alt'])] +
                                 bc.ix[i, "NOR.{}n".format(sp.ix[i, 'Alt'])])
        sp.ix[i, 'File1E'] = 0
        sp.ix[i, 'File1D'] = 0
        sp.ix[i, 'File2R'] = int(bc.ix[i, "TUM.{}p".format(sp.ix[i, 'Ref'])] +
                                 bc.ix[i, "TUM.{}n".format(sp.ix[i, 'Ref'])])
        sp.ix[i, 'File2A'] = int(bc.ix[i, "TUM.{}p".format(sp.ix[i, 'Alt'])] +
                                 bc.ix[i, "TUM.{}n".format(sp.ix[i, 'Alt'])])
        sp.ix[i, 'File2E'] = 0
        sp.ix[i, 'File2D'] = 0

    for col in sp.columns[4:]:
        sp[col] = sp[col].astype(int)

    sp.to_csv(sys.stdout, index=False)
