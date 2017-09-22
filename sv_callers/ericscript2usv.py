#!/usr/bin/env python
""" convert eriscript output to unified tsv format
"""

import sys
import argparse
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='ericscript2usv.py',
                                     description='convert ericscript output to unified sv tsv format')
    parser.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'),
                        help='ericscript output file')
    args = parser.parse_args()

    df = pd.read_table(args.input)
    df = df.rename(columns={'GeneName1': 'GeneSymbol1', 'GeneName2': 'GeneSymbol2',
                            'strand1': 'Strand1', 'strand2': 'Strand2',
                            'chr1': 'Chr1', 'chr2': 'Chr2',
                            'crossingreads': 'NumSplitReads',
                            'spanningreads': 'NumSpanningReads',
                            'homology': 'Homology',
                            'fusiontype': 'FusionType'})
    info_fields = ['Homology', 'FusionType', 'JunctionSequence', 'GeneExpr1', 'GeneExpr2', 'GeneExpr_Fused',
                   'ES', 'GJS', 'US', 'EricScore']
    df['Info'] = ''
    for i, row in df.iterrows():
        fields = []
        for field in info_fields:
            if row[field] != '' and pd.notnull(row[field]):
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
    usv_df.replace('Unable to predict breakpoint position', np.nan, inplace=True)
    usv_df.to_csv(sys.stdout, index=False, sep='\t')
