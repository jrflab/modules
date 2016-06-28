#!/usr/bin/env python2.7
""" merge base count files for facets
"""

import sys
import pandas as pd
import argparse


TUM_BC_HEADER = """
Chrom Pos Ref Alt
Refidx
TUM.TOTAL_depth
TUM.MAPQ_depth
TUM.DP
TUM.Ap TUM.Cp TUM.Gp TUM.Tp
TUM.An TUM.Cn TUM.Gn TUM.Tn
TUM.INS TUM.DEL ID
""".strip().split()

NOR_BC_HEADER = """
Chrom Pos Ref Alt
Refidx
NOR.TOTAL_depth
NOR.MAPQ_depth
NOR.DP
NOR.Ap NOR.Cp NOR.Gp NOR.Tp
NOR.An NOR.Cn NOR.Gn NOR.Tn
NOR.INS NOR.DEL ID
""".strip().split()

MERGED_HEADER = """
Chrom Pos Ref Alt
TUM.DP
TUM.Ap TUM.Cp TUM.Gp TUM.Tp
TUM.An TUM.Cn TUM.Gn TUM.Tn
NOR.DP
NOR.Ap NOR.Cp NOR.Gp NOR.Tp
NOR.An NOR.Cn NOR.Gn NOR.Tn
""".strip().split()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='facets_merge_tn.py',
                                     description='merge the facets tumor and normal base count files')
    parser.add_argument('tumor')
    parser.add_argument('normal')
    parser.add_argument('--min_cov_normal', default=25, help='minimum coverage of the normal')
    args = parser.parse_args()

    tumor_bc = pd.read_table(args.tumor, dtype={'Chrom': 'str'}, names=TUM_BC_HEADER, header=0, skiprows=0)
    normal_bc = pd.read_table(args.normal, dtype={'Chrom': 'str'}, names=NOR_BC_HEADER, header=0, skiprows=0)

    tumor_bc = tumor_bc.loc[:, tumor_bc.columns.isin(MERGED_HEADER)]
    normal_bc = normal_bc.loc[:, normal_bc.columns.isin(MERGED_HEADER)]

    normal_bc = normal_bc[(normal_bc.Ref.str.len() == 1) & (normal_bc.Alt.str.len() == 1) & (normal_bc.Ref != "N") &
                          (normal_bc['NOR.DP'] >= args.min_cov_normal)]

    tumor_normal_bc = normal_bc.merge(tumor_bc, how='inner', on=['Chrom', 'Pos', 'Ref', 'Alt'])
    tumor_normal_bc = tumor_normal_bc[MERGED_HEADER]
    #tumor_normal_bc = tumor_normal_bc.sort_values(['Chrom', 'Pos'])

    tumor_normal_bc.to_csv(sys.stdout, sep='\t', index=False)
