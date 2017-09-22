#!/usr/bin/env python
""" convert eriscript output to unified tsv format
"""

import sys
import argparse
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='starfusion2usv.py',
                                     description='convert starfusion output to unified sv tsv format')
    parser.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'),
                        help='starfusion output file')
    args = parser.parse_args()

    df = pd.read_table(args.input)
    df = df.rename(columns={'#FusionName': 'FusionName', 'JunctionReadCount': 'NumSplitReads',
                            'SpanningFragCount': 'NumSpanningReads', 'LeftGene': 'Gene1',
                            'RightGene': 'Gene2'})
    df['GeneSymbol1'], df['EnsemblGene1'] = df['Gene1'].str.split('^', 1).str
    df['GeneSymbol2'], df['EnsemblGene2'] = df['Gene2'].str.split('^', 1).str
    df['Chr1'], df['Breakpoint1'], df['Strand1'] = df['LeftBreakpoint'].str.split(':', 2).str
    df['Chr2'], df['Breakpoint2'], df['Strand2'] = df['RightBreakpoint'].str.split(':', 2).str
    df['Chr1'] = df['Chr1'].str.replace('chr', '')
    df['Chr2'] = df['Chr2'].str.replace('chr', '')
    info_fields = ['SpliceType', 'LeftBreakDinuc', 'LeftBreakEntropy', 'RightBreakDinuc', 'RightBreakEntropy','LargeAnchorSupport']
    df['Info'] = ''
    for i, row in df.iterrows():
        fields = []
        for field in info_fields:
            x = "{}=".format(field)
            if ' ' in str(row[field]):
                x += '"{}"'.format(row[field])
            else:
                x += '{}'.format(row[field])
            fields.append(x)
        df.ix[i, 'Info'] = ";".join(fields)
    usv_cols = ['Chr1', 'Breakpoint1', 'Strand1', 'Chr2', 'Breakpoint2', 'Strand2',
                'NumSplitReads', 'NumSpanningReads',
                'GeneSymbol1', 'GeneSymbol2',
                'EnsemblGene1', 'EnsemblGene2', 'Info']
    usv_df = df[usv_cols].copy()
    usv_df.to_csv(sys.stdout, index=False, sep='\t')
