#!/usr/bin/env python
"""
Plot recurrent mutations using output from jrflab modules pipeline. Outputs the
following files to given directory:

    - nr_mutations_per_gene.tsv
    - recurrent_mutations_per_gene.pdf
    - recurrent_mutations.tsv

depending on columns available in given input may also output:
    - recurrent_mutations_ccf.pdf (cancer_cell_frac column)
    - recurrent_mutations_maf.pdf (TUMOR_MAF column)
"""
import argparse
import pandas as pd
import os
import errno
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerPatch
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


class LOHLegend(mpatches.Rectangle):
    pass


class HandlerLOHLegend(HandlerPatch):
    """Handler for drawing strike through rectangle in legend"""
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        p = mpatches.Rectangle(xy=(xdescent, ydescent), width=width,
                             height=height)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        lp = mlines.Line2D((xdescent, xdescent+width), (ydescent, ydescent+height), color="k", linewidth=p._linewidth)
        lp.set_transform(trans)
        return [p, lp]


def output_recurrent_mutations(snv_fp, indel_fp, outdir):
    # parse mutations
    snv = pd.read_csv(snv_fp, dtype={"CHROM": str}, sep="\t")
    indel = pd.read_csv(indel_fp, dtype={"CHROM": str}, sep="\t")
    intersect_cols = list(set(snv.columns).intersection(indel.columns))
    muts = pd.concat([snv[intersect_cols],
                      indel[intersect_cols]],
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
        f, ax = plt.subplots()
        heatmap = sns.heatmap(gene_count.astype(bool)[gene_count.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2].T)
        f.savefig(outdir + "/recurrent_mutations_per_gene.pdf")

    # output recurrent mutations
    rec_genes = list(gene_count.astype(bool)[gene_count.astype(bool).apply(lambda x: x.sum(), axis=1) > 1].index)
    mut_split[mut_split["ANN[*].GENE_SPLIT"].isin(rec_genes)].ix[:, ["CHROM", "POS", "REF", "ALT", "TUMOR_SAMPLE", "ANN[*].GENE_SPLIT", "ANN[*].HGVS_P_SPLIT", "ANN[*].IMPACT_SPLIT", "TUMOR.DP", "TUMOR.AD"]]\
        .to_csv(outdir + "/recurrent_mutations.tsv", sep="\t", index=False)

    # recurrent identical mutations
    mut_count = muts\
        .groupby(["TUMOR_SAMPLE", "CHROM", "POS", "REF", "ALT"])["TUMOR.DP"].count().unstack("TUMOR_SAMPLE").fillna(0)
    mut_count.astype(int).to_csv(outdir + "/nr_identical_mutations.tsv", sep="\t")

    # output identical recurrent mutation plot using MAF
    if "TUMOR_MAF" in muts.columns:
        mut_maf = muts\
            .groupby(["TUMOR_SAMPLE", "ANN[*].GENE", "ANN[*].HGVS_P"])["TUMOR_MAF"].max().unstack("TUMOR_SAMPLE").fillna(0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f, ax = plt.subplots()
            heatmap = sns.heatmap(mut_maf[mut_maf.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2].T,
                mask=(mut_maf[mut_maf.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2] == 0).T,
                annot=True, vmin=0, vmax=1)
            f.savefig(outdir + "/recurrent_mutations_maf.pdf")

    # output identical recurrent mutation plot using ccf
    if "cancer_cell_frac" in muts.columns:
        mut_ccf = muts\
            .groupby(["TUMOR_SAMPLE", "ANN[*].GENE", "ANN[*].HGVS_P"])["cancer_cell_frac"].max().unstack("TUMOR_SAMPLE").fillna(0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f, ax = plt.subplots()
            heatmap = sns.heatmap(mut_ccf[mut_ccf.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2].T,
                annot=True,
                mask=(mut_ccf[mut_ccf.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2] == 0).T, vmin=0, vmax=1)
            f.savefig(outdir + "/recurrent_mutations_ccf.pdf")

    # output nonsilent mutations plot using ccf
    if all([c in muts.columns for c in "cancer_gene pathogenicity LOH clonality".split()]):
        mut_ccf = mut_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                            (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
            .drop_duplicates(subset="TUMOR_SAMPLE CHROM POS REF ALT".split())\
            .groupby(["TUMOR_SAMPLE", "ANN[*].GENE_SPLIT", "ANN[*].HGVS_P_SPLIT"])["cancer_cell_frac"].max().unstack("TUMOR_SAMPLE").fillna(0)

        # mutation annotations
        annotations = mut_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                            (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
            .groupby(["ANN[*].GENE_SPLIT", "ANN[*].HGVS_P_SPLIT"])["cancer_gene pathogenicity LOH clonality".split()].first()

        sns.set_context(context="poster", font_scale=0.2)
        sns.set_style("whitegrid")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            for sample in mut_split.TUMOR_SAMPLE.unique():
                sample_muts = mut_ccf[sample].sort_values(ascending=False)
                sample_muts = pd.DataFrame(sample_muts[sample_muts > 0])
                font_scale = min([1, 30.0/len(sample_muts)])
                sns.set_context(context="poster", font_scale=font_scale)
                f, ax = plt.subplots()
                cbar_ax = f.add_axes([0.7, .53, .15, .02])
                legend_ax = f.add_axes([0.7, .55, .15, .10])
                legend_ax.axis('off')
                heatmap = sns.heatmap(sample_muts.T, annot=False, vmin=0, vmax=1, ax=ax, cbar_ax=cbar_ax,
                                      cbar_kws={"label": "", "orientation":"horizontal"}, linewidths=.5,
                                      square=True)
                # heatmap.ax_heatmap.xaxis.set_tick_params(labelsize=3)
                heatmap.xaxis.set_label_text("")
                heatmap.xaxis.set_label_position('top')
                heatmap.xaxis.set_ticks_position('top')
                plt.setp(heatmap.xaxis.get_majorticklabels(), rotation=90)
                plt.setp(heatmap.yaxis.get_majorticklabels(), rotation=0)

                # solid cbar http://stackoverflow.com/questions/33652272
                cbar_ax.collections[0].set_edgecolor("face")

                # Add annotations
                for i, x in enumerate(heatmap.xaxis.get_ticklocs()):
                    if annotations.ix[sample_muts.index[i]].clonality == "clonal":
                        heatmap.add_patch(plt.Rectangle((x-.5, 0), 1, 1, fill=False, edgecolor='goldenrod', linewidth=2, clip_on=False))
                    if annotations.ix[sample_muts.index[i]].pathogenicity == "pathogenic":
                        heatmap.add_artist(plt.Circle((x, -1), 0.1, color='r', clip_on=False))
                    if annotations.ix[sample_muts.index[i]].cancer_gene == "true":
                        heatmap.add_artist(plt.Circle((x, -2), 0.1, color='k', clip_on=False))
                    if annotations.ix[sample_muts.index[i]].LOH == "true":
                        heatmap.add_artist(plt.Line2D((i, i + 1), (0, 1), 1, color='w'))

                # add legend
                # make room
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
                pa = plt.Line2D(range(1), range(1), markersize=5, color="w", marker='o', markerfacecolor="r")
                cg = plt.Line2D(range(1), range(1), markersize=5, color="w", marker='o', markerfacecolor="k")
                loh = LOHLegend((1, 1), 1, 1, facecolor="w", edgecolor="k", linewidth=0.1)
                clonal = mpatches.Rectangle((1, 1), 1, 1, facecolor="w", edgecolor="goldenrod", linewidth=2)
                ccf = mpatches.Patch(facecolor="w", edgecolor="w")
                legend = legend_ax.legend([cg, pa, loh, clonal, ccf],
                    ["Cancer gene", "Pathogenic", "LOH", "Clonal", "CCF"],
                    handlelength=0.8, loc='lower left',
                    handler_map={LOHLegend: HandlerLOHLegend()})
                for t in legend.get_texts():
                    t.set_ha("left")

                f.savefig(outdir + "/" + sample + "_nonsilent_mutations_ccf.pdf")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mutect_nonsynonymous", type=str, help="TSV")
    parser.add_argument("strelka_varscan_nonsynonymous", type=str, help="TSV")
    parser.add_argument("outdir", type=str, help="output directory")
    args = parser.parse_args()
    output_recurrent_mutations(args.mutect_nonsynonymous,
                                        args.strelka_varscan_nonsynonymous,
                                        args.outdir)


if __name__ == "__main__":
    main()
