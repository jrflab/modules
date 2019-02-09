#!/usr/bin/env python
""" convert integrate output to unified tsv format
"""

import sys
import argparse
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='integrate2usv.py',
                                     description='convert INTEGRATE output to unified sv tsv format')
    parser.add_argument('--breakpoints_file', type=argparse.FileType('r'),
                        help='integrate breakpoints file')
    parser.add_argument('--sum_file', type=argparse.FileType('r'),
                        help='integrate sum file')
    parser.add_argument('--exons_file', type=argparse.FileType('r'),
                        help='integrate exons file')
    args = parser.parse_args()

    df_bp = pd.read_table(args.breakpoints_file)
    df_sum = pd.read_table(args.sum_file)
    df_exon = pd.read_table(args.exons_file)
    df_bp_sum = pd.concat([df_bp, df_sum], axis=1)
    df_bp_sum = df_bp_sum.drop(['5_Prime', '3_Prime'])
    df_exon = df_exon.drop(['5P', '3P'])
    df = pd.merge(df_bp_sum, df_exon, how='left', left_on='Fusion_Candidate', right_on='#Id')

    df = df.rename(columns={'#5P': 'GeneSymbol1', '3P': 'GeneSymbol2',
                            'strand1': 'Strand1', 'strand2': 'Strand2',
                            'chr1': 'Chr1', 'chr2': 'Chr2',
                            'crossingreads': 'NumSplitReads',
                            'spanningreads': 'NumSpanningReads',
                            'homology': 'Homology',
                            'fusiontype': 'FusionType'})
    info_fields = ['Homology', 'FusionType', 'JunctionSequence', 'GeneExpr1', 'GeneExpr2', 'GeneExpr_Fused', 'ES', 'GJS', 'US', 'EricScore']
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
