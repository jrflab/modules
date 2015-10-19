#!/usr/bin/env python
"""
Plot recurrent mutations using output from jrflab modules pipeline. Outputs the
following files to given directory:

    - nr_mutations_per_gene.tsv
    - recurrent_mutations_per_gene.pdf
    - recurrent_mutations.tsv
"""
import argparse
import pandas as pd
import os
import errno
import warnings

import matplotlib
matplotlib.use('Agg')
import seaborn as sns


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def simplify(df):
    df = df.rename(columns={"#CHROM": "CHROM"})
    return df["TUMOR_SAMPLE CHROM POS ID REF ALT ANN[*].GENE ANN[*].HGVS_P ANN[*].EFFECT ANN[*].IMPACT TUMOR.DP TUMOR.AD dbNSFP_ExAC_Adj_AF dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AF dbNSFP_1000Gp3_AF G5 G5A GMAF".split()]


def split_multi_value_columns_to_records(df, columns, separator):
    assert(len(columns) > 1)
    split = pd.DataFrame(df[columns[0]].str.split(separator).tolist(), index=df.index).stack()
    split.name = "SPLIT"
    split.index = split.index.droplevel(-1)
    df_split = df.join(pd.DataFrame({c + "_SPLIT": [s for l in
                                                    df[c].str.split(separator).tolist()
                                                    for s in l] for c in
                                     columns}, index=split.index))
    df_split.index = list(range(len(df_split)))
    return df_split


def output_recurrent_mutations_per_gene(snv_fp, indel_fp, outdir):
    # parse mutations
    muts = pd.concat([simplify(pd.read_csv(snv_fp, dtype={"CHROM": str},
                                           sep="\t")),
                      simplify(pd.read_csv(indel_fp, dtype={"CHROM": str},
                                           sep="\t"))],
                     ignore_index=True)
    # reshape df to record format (columns contain | separated values)
    mut_split = split_multi_value_columns_to_records(muts,
                                                    "ANN[*].GENE ANN[*].HGVS_P ANN[*].EFFECT ANN[*].IMPACT".split(),
                                                    "|")
    # create output dir
    mkdir_p(outdir)

    # mutation per gene count
    gene_count = mut_split[mut_split["ANN[*].IMPACT_SPLIT"].isin(["MODERATE", "HIGH"])]\
        .groupby(["TUMOR_SAMPLE", "ANN[*].GENE_SPLIT"]).CHROM.count().unstack("TUMOR_SAMPLE").fillna(0)
    gene_count.astype(int).to_csv(outdir + "/nr_mutations_per_gene.tsv", sep="\t")

    # output recurrent mutation plot
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        heatmap = sns.clustermap(gene_count.astype(bool)[gene_count.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2],
            col_cluster=False, row_cluster=False)
        heatmap.savefig(outdir + "/recurrent_mutations_per_gene.pdf")

    # output recurrent mutations
    rec_genes = list(gene_count.astype(bool)[gene_count.astype(bool).apply(lambda x: x.sum(), axis=1) > 1].index)
    mut_split[mut_split["ANN[*].GENE_SPLIT"].isin(rec_genes)].ix[:, ["CHROM", "POS", "REF", "ALT", "TUMOR_SAMPLE", "ANN[*].GENE_SPLIT", "ANN[*].HGVS_P_SPLIT", "ANN[*].IMPACT_SPLIT", "TUMOR.DP", "TUMOR.AD"]]\
        .to_csv(outdir + "/recurrent_mutations.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mutect_nonsynonymous", type=str, help="TSV")
    parser.add_argument("strelka_varscan_nonsynonymous", type=str, help="TSV")
    parser.add_argument("outdir", type=str, help="output directory")
    args = parser.parse_args()
    output_recurrent_mutations_per_gene(args.mutect_nonsynonymous,
                                        args.strelka_varscan_nonsynonymous,
                                        args.outdir)


if __name__ == "__main__":
    main()
