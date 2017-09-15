#!/usr/bin/env python
"""
Melt a somatic vcf using vcf_melt. Split and merge the table on SAMPLE.
Create a table with one line per variant
"""
import argparse
import sys
import pandas as pd
import subprocess
import re


def add_maf(df):
    rv = df.copy()

    def f(x):
        if "," in str(x["AD.TUMOR"]) and sum(map(float, x["AD.TUMOR"].split(','))) > 0:
            return float(x["AD.TUMOR"].split(",")[1]) / sum(map(float, x["AD.TUMOR"].split(',')))
        else:
            return float('nan')

    def g(x):
        if "," in str(x["AD.NORMAL"]) and sum(map(float, x["AD.NORMAL"].split(','))) > 0:
            return float(x["AD.NORMAL"].split(",")[1]) / sum(map(float, x["AD.NORMAL"].split(',')))
        else:
            return float('nan')

    if len(df) > 0 and all([x in df.columns for x in ['AD.TUMOR', 'AD.NORMAL']]):
        rv["TUMOR_MAF"] = df.apply(f, axis=1)
        rv["NORMAL_MAF"] = df.apply(g, axis=1)
    elif len(df) > 0 and all([x in df.columns for x in ['NV.TUMOR', 'NV.NORMAL', 'NR.TUMOR', 'NR.NORMAL']]):
        rv["TUMOR_MAF"] = df['NV.TUMOR'].apply(lambda x: int(str(x).split(',')[0])) / \
            df['NR.TUMOR'].apply(lambda x: int(str(x).split(',')[0]))
        rv["NORMAL_MAF"] = df['NV.NORMAL'].apply(lambda x: int(str(x).split(',')[0])) / \
            df['NR.NORMAL'].apply(lambda x: int(str(x).split(',')[0]))
    else:
        rv["TUMOR_MAF"] = pd.Series()
        rv["NORMAL_MAF"] = pd.Series()
    return rv


def add_dp(df):
    rv = df.copy()

    def h(x):
        if "," in str(x["AD.TUMOR"]):
            return sum(map(float, x["AD.TUMOR"].split(',')))
        else:
            return float('nan')

    def i(x):
        if "," in str(x["AD.NORMAL"]):
            return sum(map(float, x["AD.NORMAL"].split(',')))
        else:
            return float('nan')

    if len(df) > 0 and all([x in df.columns for x in ['AD.TUMOR', 'AD.NORMAL']]):
        rv["TUMOR_DP"] = df.apply(h, axis=1)
        rv["NORMAL_DP"] = df.apply(i, axis=1)
    elif len(df) > 0 and all([x in df.columns for x in ['NR.TUMOR', 'NR.NORMAL']]):
        rv["TUMOR_DP"] = df['NR.TUMOR'].apply(lambda x: int(str(x).split(',')[0]))
        rv["NORMAL_DP"] = df['NR.NORMAL'].apply(lambda x: int(str(x).split(',')[0]))
    else:
        rv["TUMOR_DP"] = pd.Series()
        rv["NORMAL_DP"] = pd.Series()
    return rv


def merge_ann_cols(df):
    rv = df.copy()
    if len(df) > 0:
        if 'MutationTaster_pred' in rv and 'MT_pred' in rv:
            rv.ix[rv['is_indel'] & rv['MutationTaster_pred'].isnull(), 'MutationTaster_pred'] = \
                rv.ix[rv['is_indel'] & rv['MutationTaster_pred'].isnull(), 'MT_pred']
            rv.ix[rv['MutationTaster_pred'].astype('str').str.contains('disease').fillna(False), 'MutationTaster_pred'] = 'D'
            rv.ix[rv['MutationTaster_pred'].astype('str').str.contains('poly').fillna(False), 'MutationTaster_pred'] = 'N'
        if 'dbNSFP_PROVEAN_pred' in rv:
            rv.ix[rv['dbNSFP_PROVEAN_pred'].astype('str').str.contains('N').fillna(False), 'PROVEAN_pred'] = 'N'
            rv.ix[rv['dbNSFP_PROVEAN_pred'].astype('str').str.contains('D').fillna(False), 'PROVEAN_pred'] = 'D'
        if 'provean_pred' in rv:
            rv.ix[rv['is_indel'] & rv['provean_pred'].astype('str').str.contains('N').fillna(False),
                  'PROVEAN_pred'] = 'N'
            rv.ix[rv['is_indel'] & rv['provean_pred'].astype('str').str.contains('D').fillna(False),
                  'PROVEAN_pred'] = 'D'
    return rv


def add_cancer_gene(df):
    rv = df.copy()
    if len(df) > 0 and all([x in df.columns for x in ['cancer_gene_census', 'kandoth', 'lawrence']]):
        rv["cancer_gene"] = rv['cancer_gene_census'] | rv['kandoth'] | rv['lawrence']
        rv["num_cancer_gene"] = df[['cancer_gene_census', 'kandoth', 'lawrence']].sum(axis=1)
    return rv

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--normal', required=True, help='Normal sample')
    parser.add_argument('vcf_file')

    args = parser.parse_args()

    if args.vcf_file.endswith('gz'):
        melt = subprocess.Popen('zcat {} | vcf_melt'.format(args.vcf_file), stdout=subprocess.PIPE, shell=True)
    else:
        melt = subprocess.Popen('vcf_melt < {}'.format(args.vcf_file), stdout=subprocess.PIPE, shell=True)
    df = pd.read_table(melt.stdout, sep='\t', low_memory=False, na_values=['None', '.'])
    if 'info.REF' in df:
        df = df.drop('info.REF', axis=1)
    df.ALT = df.ALT.str.replace(r'[\[\]]', r'')
    df.ALT = df.ALT.str.replace(r',.*', r'')
    df['is_indel'] = (df.ALT.str.len() > 1) | (df.REF.str.len() > 1)
    tdf = df.loc[df.SAMPLE != args.normal]
    ndf = df.loc[df.SAMPLE == args.normal]

    merge_cols = ['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'is_indel']
    info_cols = tdf.columns[pd.Series(tdf.columns).astype('str').str.startswith('info.')]
    merge_cols.extend(info_cols.tolist())
    mdf = pd.merge(tdf, ndf, how='left', on=merge_cols, suffixes=('.TUMOR', '.NORMAL'))
    mdf.rename(columns=lambda x: re.sub('info.', '', x), inplace=True)

    mdf = add_maf(mdf)
    mdf = add_dp(mdf)
    mdf = merge_ann_cols(mdf)
    mdf = add_cancer_gene(mdf)

    chasm_pred_columns = [c for c in mdf.columns if "chasm_pred" in c]
    summary_columns = ("CHROM,POS,REF,ALT,SAMPLE.TUMOR,SAMPLE.NORMAL,"
                       "variantCaller,UPS_Coord,SYMBOL,Variant_Classification,HGVSp_Short,"
                       "TUMOR_MAF,NORMAL_MAF,TUMOR_DP,NORMAL_DP,offTarget,"
                       "fuentes,dgd,oncoKB_level,oncoKB_cancer_type,"
                       "cancer_gene_census,kandoth,lawrence,num_cancer_gene,hap_insuf,"
                       "AF,ExAC_AF,MutationTaster_pred,PROVEAN_pred,FATHMM_pred," +
                       ",".join(chasm_pred_columns) + ',' +
                       "facetsLOHCall,parssnp_pred,pathogenicity,HOTSPOT,HOTSPOT_INTERNAL,cmo_hotspot," +
                       "clonalStatus,ccf,samples_called_in").split(",")
    cols = [x for x in summary_columns if x in mdf]

    mdf.ix[:, cols].to_csv(sys.stdout, sep='\t', index=False)
