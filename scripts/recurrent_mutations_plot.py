#!/usr/bin/env python
"""
Plot recurrent mutations using output from jrflab modules pipeline. Outputs the
following files to given directory:

    - nr_mutations_per_gene.tsv
    - recurrent_mutations_per_gene.pdf
    - recurrent_mutations.tsv

depending on columns available in given input may also output:
    - recurrent_mutations_ccf.pdf (ccf column)
    - recurrent_mutations_maf.pdf (TUMOR_MAF column)
"""
import argparse
import pandas as pd
import numpy as np
import os
import errno
import warnings
import yaml
import re
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerPatch
import seaborn as sns


def three_to_one_amino_acid_code(x):
    # 3 > 1 AA code mapping
    d = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
         'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
         'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
         'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}
    for aa in d:
        x = re.sub(aa, d[aa], x)
    return x


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
    return df["TUMOR_SAMPLE CHROM POS ID REF ALT ANN[*].GENE ANN[*].HGVS_P ANN[*].EFFECT ANN[*].IMPACT TUMOR_DP TUMOR.AD dbNSFP_ExAC_Adj_AF dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AF dbNSFP_1000Gp3_AF G5 G5A GMAF".split()]


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


def make_most_severe_effect_column(mut_split):
    """Normalize ANN[*].EFFECT_SPLIT column to most severe mutation effect. If
    the mutation is a hotspot that is even more severe. The MOST_SEVERE_EFFECT
    column is returned as a column of integers, where the highest number is the
    most severe. The mapping from severity number to effect name is also returned"""
    rv = mut_split.copy()
    # make ordered dict of normalization patterns so when normalization is done
    # in order of keys, the most severe change is applied last
    normalize = OrderedDict([
        ("Silent", ".*(silent|synonymous|intron|intragenic|utr|igr|rna|non_coding|noncoding).*"),
        ("Up/downstream variant", ".*(upstream|downstream).*"),
        ("Start or Stop codon change", ".*(start|stop).*"),
        ("Splice site variant", ".*splice.*"),
        ("Missense SNV",  ".*(missense|non_synonymous|nonsynonymous).*"),
        ("Inframe In-Del",  ".*(codon_deletion|codon_insertion|inframe|in_frame).*"),
        ("Frameshift In-Del", ".*shift.*"),
        ("Truncating SNV", ".*(trun|stop).*"),
    ])
    rv["MOST_SEVERE_EFFECT"] = [np.nan] * len(rv["ANN[*].EFFECT_SPLIT"])

    for i, (effect, pattern) in enumerate(normalize.items()):
        rv.loc[rv["ANN[*].EFFECT_SPLIT"].str.match(re.compile(pattern, re.IGNORECASE)).astype(bool), 'MOST_SEVERE_EFFECT'] = i
    rv.loc[rv["HOTSPOT"] == True, 'MOST_SEVERE_EFFECT'] = len(normalize)

    number_effect_mapping = dict(enumerate(normalize.keys()))
    number_effect_mapping[len(normalize)] = "Hotspot"

    return rv, number_effect_mapping


def output_recurrent_mutations_type_single_figure_plot_only(mut_eff, gene_order, sample_order, number_effect_mapping, outputfile):
    # create color map
    colors = list(reversed(sns.color_palette("Set2", len(number_effect_mapping) - 1)))
    colors += [(1.0, 0.0, 0.0)]
    cm = matplotlib.colors.ListedColormap(colors)

    # draw the plot
    sns.set_context("notebook")
    width = len(mut_eff.columns)
    height = len(mut_eff)*.5
    f, ax = plt.subplots(figsize=(width, height))

    legend_ax = f.add_axes([0.95, .875, 1/width, 1/height])
    legend_ax.axis('off')
    heatmap = sns.heatmap(mut_eff.ix[gene_order, sample_order],
        ax=ax,
        cbar=False,
        square=True,
        annot=False,
        cmap=cm,
        vmin=0,
        vmax=len(colors),
        mask=(mut_eff.ix[gene_order, sample_order] < 0),
        linewidths=.5)
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_label_text("")
    ax.xaxis.tick_top()
    ax.yaxis.set_label_text("")
    plt.setp(heatmap.xaxis.get_majorticklabels(), rotation=90)

    # add color map to legend
    patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]
    legend = legend_ax.legend(patches,
        [number_effect_mapping[k] for k in sorted(number_effect_mapping.keys())],
        handlelength=0.8, loc='lower left')
    for t in legend.get_texts():
        t.set_ha("left")

    f.savefig(outputfile, bbox_inches='tight')


def output_recurrent_mutations_type_single_figure(mut_split, gene_order, sample_order, outputfile, per_mutation=False):
    if per_mutation:
        groupby = "ANN[*].GENE_SPLIT ANN[*].HGVS_P_SPLIT_SHORT".split()
    else:
        groupby = "ANN[*].GENE_SPLIT".split()  # plot per gene
    mut_eff_split, number_effect_mapping = make_most_severe_effect_column(mut_split)
    mut_eff = mut_eff_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                                                  (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
            .groupby(["TUMOR_SAMPLE"] + groupby)["MOST_SEVERE_EFFECT"].max().unstack("TUMOR_SAMPLE").fillna(-1).astype(int)
    output_recurrent_mutations_type_single_figure_plot_only(mut_eff, gene_order, sample_order, number_effect_mapping, outputfile)


def output_recurrent_mutations_ccf_per_mutation_single_figure(mut_split, gene_order, sample_order, outputfile, plot_maf=False):
    output_recurrent_mutations_ccf_single_figure_plot_only(mut_split, gene_order, sample_order, outputfile, plot_maf, per_mutation=True)


def output_recurrent_mutations_ccf_single_figure_plot_only(mut_split, gene_order, sample_order, outputfile, plot_maf=False, per_mutation=False):
    # ccf per gene
    ccf_column = "TUMOR_MAF" if plot_maf else "ccf"
    mut_ccf, mut_annotations, mut_ann_per_sample = get_mutation_heatmap_dfs(mut_split, ccf_column, per_mutation=per_mutation)

    # draw the plot
    sns.set_context("notebook")
    width = len(mut_ccf.columns)
    height = len(mut_ccf)*.5
    f, ax = plt.subplots(figsize=(width, height))

    cbar_ax = f.add_axes([0.97, .874, 1.5/width, 0.10/height])
    legend_ax = f.add_axes([0.95, .875, 1/width, 1/height])
    legend_ax.axis('off')

    cmap = None if plot_maf else "Blues"

    heatmap = sns.heatmap(mut_ccf.ix[gene_order, sample_order],
        ax=ax,
        cbar_ax=cbar_ax,
        cbar_kws={"label": "", "orientation": "horizontal"},
        square=True,
        annot=False,
        cmap=cmap,
        mask=(mut_ccf.ix[gene_order, sample_order] == 0),
        linewidths=.5,
        vmin=0,
        vmax=1)
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_label_text("")
    ax.xaxis.tick_top()
    ax.yaxis.set_label_text("")
    plt.setp(heatmap.xaxis.get_majorticklabels(), rotation=90)

    # solid cbar http://stackoverflow.com/questions/33652272
    cbar_ax.collections[0].set_edgecolor("face")

    # add mutation annotations
    for i, y in enumerate(heatmap.yaxis.get_ticklocs()):
        if mut_annotations.ix[mut_ccf.ix[gene_order].index[i]].cancer_gene_census:
            heatmap.add_artist(matplotlib.patches.Ellipse((len(sample_order)+.2, len(mut_ccf.ix[gene_order])-y), .2,
                                                        .2,
                                                        color='k', clip_on=False))
        if mut_annotations.ix[mut_ccf.ix[gene_order].index[i]].pathogenic:
            heatmap.add_artist(matplotlib.patches.Ellipse((len(sample_order)+.2+.5, len(mut_ccf.ix[gene_order])-y), .2,
                                                        .2,
                                                        color='r', clip_on=False))
    # add per sample mutation annotations
    for i, x in enumerate(heatmap.xaxis.get_ticklocs()):
        for j, y in enumerate(heatmap.yaxis.get_ticklocs()):
            try:
                ann = mut_ann_per_sample.ix[mut_ccf.ix[gene_order].index[j]].ix[sample_order[i]]
            except KeyError:
                continue
            if not plot_maf and ann.clonal:
                heatmap.add_patch(
                    plt.Rectangle((x-.5, len(mut_ccf.ix[gene_order])-y-.5),
                                1, 1, fill=False, edgecolor='goldenrod', linewidth=2, clip_on=False))
            if not plot_maf and ann.facetsLOHCall:
                heatmap.add_artist(plt.Line2D((i, i + 1), (len(mut_ccf.ix[gene_order])-j-1, len(mut_ccf.ix[gene_order])-j), 1, color='w'))

    # add legend
    # box = ax.get_position() make room
    # ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    cg = plt.Line2D(list(range(1)), list(range(1)), markersize=5, color="w", marker='o', markerfacecolor="k")
    pa = plt.Line2D(list(range(1)), list(range(1)), markersize=5, color="w", marker='o', markerfacecolor="r")
    loh = LOHLegend((1, 1), 1, 1, facecolor="w", edgecolor="k", linewidth=1)
    clonal = mpatches.Rectangle((1, 1), 1, 1, facecolor="w", edgecolor="goldenrod", linewidth=2)
    ccf = mpatches.Patch(facecolor="w", edgecolor="w")
    ccf_text = "MAF" if plot_maf else "CCF"
    legend = legend_ax.legend([cg, pa, loh, clonal, ccf],
        ["Cancer gene", "Pathogenic", "LOH", "Clonal", ccf_text],
        handlelength=1, loc='lower left', fontsize=13,
        handler_map={LOHLegend: HandlerLOHLegend()})
    for t in legend.get_texts():
        t.set_ha("left")

    f.savefig(outputfile, bbox_inches='tight')


def get_mutation_heatmap_dfs(mut_split, ccf_column, per_mutation=False):
    if per_mutation:
        groupby = "ANN[*].GENE_SPLIT ANN[*].HGVS_P_SPLIT_SHORT".split()
    else:
        groupby = "ANN[*].GENE_SPLIT".split()  # plot per gene
    mut_ccf = mut_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                                                  (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
            .drop_duplicates(subset="TUMOR_SAMPLE CHROM POS REF ALT".split())\
            .groupby(["TUMOR_SAMPLE"] + groupby)[ccf_column].max().unstack("TUMOR_SAMPLE").fillna(0)
    # mutation annotations
    mut_annotations = mut_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                                                     (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
            .groupby(groupby)["cancer_gene_census facetsLOHCall pathogenicity clonalStatus".split()].apply(
                    lambda x: pd.Series({"cancer_gene_census": len(x[x.cancer_gene_census == True]) > 0,
                                         "pathogenic": ((x.pathogenicity == 'pathogenic') | (x.pathogenicity == 'likely_pathogenic')).sum() > 0}))
    # gene annotations per mutation per sample
    mut_ann_per_sample = mut_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                        (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
        .groupby(groupby + ["TUMOR_SAMPLE"])["facetsLOHCall clonalStatus".split()].apply(
            lambda x: pd.Series({"facetsLOHCall": len(x[x.facetsLOHCall == True]) > 0,
                                "clonal": len(x[x.clonalStatus == "clonal"]) > 0,
                                }))

    return mut_ccf, mut_annotations, mut_ann_per_sample


def get_heatmap_order(gene_count, gene_annotations, mut_ccf_genes, ccf_column='ccf', config=None):
    # sort recurrently mutated genes by number of hits in each sample
    # 1. number of hits in each sample
    # 2. ccf
    # 3. cancer_gene/pathogenic
    top_recurrent = gene_count.astype(bool).apply(lambda x: x.sum(), axis=1).sort_values(ascending=False)
    top_recurrent_genes_order_df = \
        pd.DataFrame({"gene":
                      top_recurrent[top_recurrent > 1].index,
                      "top_recurrent":
                      top_recurrent[top_recurrent > 1],
                      "ccf_sum":
                      [mut_ccf_genes.loc[g].sum() for i, g in enumerate(top_recurrent[top_recurrent > 1].index)],
                      "cancer_gene_census":
                      [gene_annotations.loc[g].cancer_gene_census for i, g in enumerate(top_recurrent[top_recurrent > 1].index)],
                      "pathogenic":
                      [gene_annotations.loc[g].pathogenic for i, g in enumerate(top_recurrent[top_recurrent > 1].index)]
                    })
    top_recurrent_genes_order = list(top_recurrent_genes_order_df.sort_values(["top_recurrent", "ccf_sum", "cancer_gene_census", "pathogenic", "gene"], ascending=[False, False, False, False, True]).gene)
    assert(len(list(top_recurrent[top_recurrent > 1].index)) == len(top_recurrent_genes_order))

    if config and 'sample_order' in config:
        sample_order = config['sample_order']
    else:
        # order sample by
        # 1. having a hit in most recurrently mutated gene
        # 2. nr of mutations in recurrently mutated genes
        # dataframe cols: sample top_recurrent nr_hits_recurrent
        sample_order_df = \
            pd.DataFrame({"sample":
                            [s for s in gene_count.columns],
                            "top_recurrent":
                            [min([top_recurrent.index.get_loc(g) for g in list(gene_count[gene_count.loc[:,s] > 0].index)]) for s in gene_count.columns],
                            "nr_hits_recurrent":
                            [gene_count.ix[top_recurrent[top_recurrent > 1].index, s].astype(bool).sum() for s in gene_count.columns]
            })
        sample_order = list(sample_order_df.sort_values(["top_recurrent", "nr_hits_recurrent", "sample"],
                                                        ascending=[True, False, True])["sample"])

    # order gene by
    # 1. most recurrent
    # 2. sample column order
    # 3. ccf
    sample_order_dict = {v:i for i, v in enumerate(sample_order)}

    # determine non recurrent genes order
    # gene top_sample_column ccf

    # get the min sample column index for each gene
    top_sample_column = [min([sample_order_dict[s] for s in gene_count.ix[g][gene_count.ix[g] > 0].index])
                                                   for g in top_recurrent[top_recurrent == 1].index]
    nonrecurrent_genes_order_df = \
            pd.DataFrame({"gene":
                           top_recurrent[top_recurrent == 1].index,
                          "top_sample_column":
                           top_sample_column,
                          "ccf":
                           [mut_ccf_genes.ix[g, sample_order[top_sample_column[i]]] for i, g in enumerate(top_recurrent[top_recurrent == 1].index)]
            })
    non_recurrent_genes_order = list(nonrecurrent_genes_order_df.sort_values(["top_sample_column", "ccf", "gene"], ascending=[True, False, False]).gene)
    gene_order = top_recurrent_genes_order + non_recurrent_genes_order

    return gene_order, sample_order



def output_recurrent_mutations_ccf_per_gene_single_figure(mut_split, gene_count, outputfile, plot_maf=False, config=None):
    """Plot the standard CCF plot per gene with LOH and cancer gene annotations."""
    # ccf per gene/mutation
    ccf_column = "TUMOR_MAF" if plot_maf else "ccf"
    mut_ccf_genes, gene_annotations, mut_ann_per_sample = get_mutation_heatmap_dfs(mut_split, ccf_column, per_mutation=False)
    gene_order, sample_order = get_heatmap_order(gene_count, gene_annotations, mut_ccf_genes, config=config)

    output_recurrent_mutations_ccf_single_figure_plot_only(mut_split, gene_order, sample_order, outputfile, plot_maf, per_mutation=False)

    # TODO: Should rewrite signature of this function or maybe remove entirely.
    # One can use get_mutation_heatmap_dfs/get_heatmap_order now to get the
    # order. No need to return this
    return gene_order, sample_order, gene_annotations


def output_recurrent_mutations(snv_fp, indel_fp, outdir, summary_config=None):
    if summary_config:
        config = yaml.load(open(summary_config))

    # parse mutations
    snv = pd.read_csv(snv_fp, dtype={"CHROM": str}, sep="\t")
    if 'ccf' in snv.columns:
        snv['ccf'] = snv['ccf'].replace({'.':np.nan}).astype(float)
    indel = pd.read_csv(indel_fp, dtype={"CHROM": str}, sep="\t")
    if 'ccf' in indel.columns:
        indel['ccf'] = indel['ccf'].replace({'.':np.nan}).astype(float)
    intersect_cols = list(set(snv.columns).intersection(indel.columns))
    muts = pd.concat([snv[intersect_cols],
                      indel[intersect_cols]],
                     ignore_index=True)
    # reshape df to record format (columns contain | separated values)
    mut_split = split_multi_value_columns_to_records(muts,
                                                    "ANN[*].GENE ANN[*].HGVS_P ANN[*].EFFECT ANN[*].IMPACT".split(),
                                                    "|")
    # TODO: this might be better moved to mutation_summary_excel, so there is
    # an extra column with the short AA codes
    mut_split["ANN[*].HGVS_P_SPLIT_SHORT"] = mut_split["ANN[*].HGVS_P_SPLIT"].apply(lambda x: three_to_one_amino_acid_code(x[2:]) if str(x).startswith('p.') else x)

    # create output dir
    mkdir_p(outdir)

    # mutation per gene count
    gene_count = mut_split[mut_split["ANN[*].IMPACT_SPLIT"].isin(["MODERATE", "HIGH"]) &
                           (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
        .drop_duplicates(subset="TUMOR_SAMPLE CHROM POS REF ALT".split())\
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
    mut_split[mut_split["ANN[*].GENE_SPLIT"].isin(rec_genes)].ix[:, ["CHROM", "POS", "REF", "ALT", "TUMOR_SAMPLE", "ANN[*].GENE_SPLIT", "ANN[*].HGVS_P_SPLIT", "ANN[*].IMPACT_SPLIT", "TUMOR_DP", "TUMOR.AD"]]\
        .to_csv(outdir + "/recurrent_mutations.tsv", sep="\t", index=False)

    # recurrent identical mutations
    mut_count = muts\
        .groupby(["TUMOR_SAMPLE", "CHROM", "POS", "REF", "ALT"])["TUMOR_DP"].count().unstack("TUMOR_SAMPLE").fillna(0)
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
    if "ccf" in muts.columns:
        mut_ccf = muts\
            .groupby(["TUMOR_SAMPLE", "ANN[*].GENE", "ANN[*].HGVS_P"])["ccf"].max().unstack("TUMOR_SAMPLE").fillna(0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f, ax = plt.subplots()
            heatmap = sns.heatmap(mut_ccf[mut_ccf.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2].T,
                annot=True,
                mask=(mut_ccf[mut_ccf.astype(bool).apply(lambda x: x.sum(), axis=1) >= 2] == 0).T, vmin=0, vmax=1)
            f.savefig(outdir + "/recurrent_mutations_ccf.pdf")

    # output nonsilent mutations plot using ccf if all required annotations are available
    if all([c in muts.columns for c in "cancer_gene_census pathogenicity facetsLOHCall clonalStatus".split()]):

        # plot combined ccf plot for all mutations per gene with annotations
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gene_order, sample_order, gene_annotations = output_recurrent_mutations_ccf_per_gene_single_figure(mut_split, gene_count, outdir + "/nonsilent_mutations_ccf.pdf", config=config)
            pathogenic_cgc_genes = set(gene_annotations[gene_annotations.cancer_gene_census | gene_annotations.pathogenic].index)
            output_recurrent_mutations_ccf_per_gene_single_figure(mut_split[mut_split['ANN[*].GENE_SPLIT'].isin(pathogenic_cgc_genes)],
                                                                  gene_count.loc[pathogenic_cgc_genes],
                                                                  outdir + "/nonsilent_mutations_ccf_pathogenic_cgc_genes.pdf", config=config)
            output_recurrent_mutations_ccf_per_mutation_single_figure(mut_split, gene_order, sample_order, outdir + "/nonsilent_mutations_ccf_per_mutation.pdf")
            output_recurrent_mutations_ccf_per_mutation_single_figure(mut_split, gene_order, sample_order, outdir + "/nonsilent_mutations_maf_per_mutation.pdf", plot_maf=True)
            output_recurrent_mutations_ccf_per_mutation_single_figure(mut_split[mut_split['ANN[*].GENE_SPLIT'].isin(pathogenic_cgc_genes)],
                                                                      gene_order, sample_order, outdir + "/nonsilent_mutations_ccf_pathogenic_cgc_genes_per_mutation.pdf")
            output_recurrent_mutations_ccf_per_gene_single_figure(mut_split, gene_count, outdir + "/nonsilent_mutations_maf.pdf", plot_maf=True, config=config)
            output_recurrent_mutations_type_single_figure(mut_split, gene_order, sample_order, outdir + "/nonsilent_mutations_type.pdf")
            output_recurrent_mutations_type_single_figure(mut_split[mut_split['ANN[*].GENE_SPLIT'].isin(pathogenic_cgc_genes)],
                                                                   [g for g in gene_order if g in pathogenic_cgc_genes],
                                                                   sample_order,
                                                                  outdir + "/nonsilent_mutations_type_pathogenic_cgc_genes.pdf")
            output_recurrent_mutations_type_single_figure(mut_split, gene_order, sample_order, outdir + "/nonsilent_mutations_type_per_mutation.pdf", per_mutation=True)
            output_recurrent_mutations_type_single_figure(mut_split[mut_split['ANN[*].GENE_SPLIT'].isin(pathogenic_cgc_genes)],
                                                                   [g for g in gene_order if g in pathogenic_cgc_genes],
                                                                   sample_order,
                                                                  outdir + "/nonsilent_mutations_type_pathogenic_cgc_genes_per_mutation.pdf", per_mutation=True)


        mut_ccf = mut_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                            (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
            .drop_duplicates(subset="TUMOR_SAMPLE CHROM POS REF ALT".split())\
            .groupby(["TUMOR_SAMPLE", "ANN[*].GENE_SPLIT", "ANN[*].HGVS_P_SPLIT"])["ccf"].max().unstack("TUMOR_SAMPLE").fillna(0)

        # mutation annotations
        annotations = mut_split[(mut_split["ANN[*].IMPACT_SPLIT"].isin(["HIGH", "MODERATE"])) &
                            (mut_split["ANN[*].HGVS_P_SPLIT"].str.contains("p."))]\
            .groupby(["ANN[*].GENE_SPLIT", "ANN[*].HGVS_P_SPLIT"])["cancer_gene_census pathogenicity facetsLOHCall clonalStatus".split()].first()

        sns.set_context(context="poster", font_scale=0.2)
        sns.set_style("whitegrid")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            for sample in mut_split.TUMOR_SAMPLE.unique():
                sample_muts = mut_ccf[sample].sort_values(ascending=False)
                sample_muts = pd.DataFrame(sample_muts[sample_muts > 0])
                #font_scale = min([1, 30.0/len(sample_muts)])
                sns.set_context(context="notebook")
                width = len(sample_muts) * 2
                height = 10

                f, ax = plt.subplots(figsize=(width, height))
                cbar_ax = f.add_axes([0.7, .53, 1.5/width, .10/height])
                legend_ax = f.add_axes([0.7, .55, 1/width, 1/height])
                legend_ax.axis('off')
                heatmap = sns.heatmap(sample_muts.T, annot=False, vmin=0, vmax=1, ax=ax, cbar_ax=cbar_ax,
                                      cmap="Blues",
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
                    if annotations.ix[sample_muts.index[i]].clonalStatus == "clonal":
                        heatmap.add_patch(plt.Rectangle((x-.5, 0), 1, 1, fill=False, edgecolor='goldenrod', linewidth=2, clip_on=False))
                    if annotations.ix[sample_muts.index[i]].pathogenicity == "pathogenic":
                        heatmap.add_artist(plt.Circle((x, -.5), .5/height, color='r', clip_on=False))
                    if annotations.ix[sample_muts.index[i]].cancer_gene_census == True:
                        heatmap.add_artist(plt.Circle((x, -1), .5/height, color='k', clip_on=False))
                    if annotations.ix[sample_muts.index[i]].facetsLOHCall == True:
                        heatmap.add_artist(plt.Line2D((i, i + 1), (0, 1), 1, color='w'))

                # add legend
                # make room
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
                pa = plt.Line2D(list(range(1)), list(range(1)), markersize=5, color="w", marker='o', markerfacecolor="r")
                cg = plt.Line2D(list(range(1)), list(range(1)), markersize=5, color="w", marker='o', markerfacecolor="k")
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
    parser.add_argument("--summary_config", default=None, type=str,
        help="Summary config yaml for summary/plotting options")
    args = parser.parse_args()
    output_recurrent_mutations(args.mutect_nonsynonymous,
                                        args.strelka_varscan_nonsynonymous,
                                        args.outdir,
                                        args.summary_config)


if __name__ == "__main__":
    main()
