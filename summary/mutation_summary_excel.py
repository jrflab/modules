#!/usr/bin/env python
"""
Create mutation summary excel from given tsvs
"""
import argparse
import pandas as pd
import os
import errno


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


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


def filter_annotations_with_impact(df, impact, effect=".*"):
    """Get only annotations with given impact"""
    # check if df is empty
    if len(df) == 0:
        return df

    def merge_impact(muts, sep):
        return pd.Series({"ANN[*].EFFECT_MERGE": sep.join(muts["ANN[*].EFFECT_SPLIT"]),
                             "ANN[*].IMPACT_MERGE": sep.join(muts["ANN[*].IMPACT_SPLIT"]),
                             "ANN[*].HGVS_P_MERGE": sep.join(muts["ANN[*].HGVS_P_SPLIT"]),
                             "ANN[*].HGVS_C_MERGE": sep.join(muts["ANN[*].HGVS_C_SPLIT"]),
                             "ANN[*].GENE_MERGE": sep.join(muts["ANN[*].GENE_SPLIT"])})

    # assume one unique entry per variant in df
    assert(len(df.drop_duplicates("TUMOR_SAMPLE CHROM POS REF ALT".split())) == len(df))
    # reshape df to record format (columns contain | separated values)
    df_split = split_multi_value_columns_to_records(df,
                                                    "ANN[*].GENE ANN[*].HGVS_P ANN[*].HGVS_C ANN[*].EFFECT ANN[*].IMPACT".split(),
                                                    "|")
    # select only variants with given impact/effect
    df_split = df_split[df_split["ANN[*].IMPACT_SPLIT"].str.match(impact) & df_split["ANN[*].EFFECT_SPLIT"].str.match(effect)]
    # merge remaining variants back to | separated values
    imp_sel = df_split.groupby("TUMOR_SAMPLE CHROM POS REF ALT".split()).apply(lambda x: merge_impact(x, "|"))
    if imp_sel.empty:
        return imp_sel
    rv = df.set_index("TUMOR_SAMPLE CHROM POS REF ALT".split()).join(imp_sel,
        how="inner")
    for c in "ANN[*].GENE ANN[*].HGVS_P ANN[*].HGVS_C ANN[*].EFFECT ANN[*].IMPACT".split():
        rv[c] = rv[c + "_MERGE"]
        del rv[c + "_MERGE"]

    return rv.reset_index()


def create_absolute_df(absolute_somatic_txts, absolute_segments):
    absdf = pd.concat([pd.read_csv(asegs, sep="\t", dtype={"Chromosome": str}) for asegs
               in absolute_somatic_txts], ignore_index=True)
    absdf["Chromosome"] = absdf.Chromosome.replace({"23": "X", "24": "Y"})
    absdf["TUMOR_SAMPLE"] = absdf.Sample.apply(lambda x: x.split("_")[0])
    absdf["NORMAL_SAMPLE"] = absdf.Sample.apply(lambda x: x.split("_")[1])
    absdf.drop([c for c in absdf.columns if c not in "TUMOR_SAMPLE NORMAL_SAMPLE Ref Alt Chromosome Position".split()],
                         axis=1,
                         inplace=True)
    absdf.rename(columns={"Ref": "REF", "Alt": "ALT", "Chromosome": "CHROM", "Position": "POS"}, inplace=True)
    abscolumns = "cancer_cell_frac Pr_somatic_clonal ccf_CI95_low".split()
    absegdf = pd.concat([pd.read_csv(f, dtype={"Chromosome": str}, sep="\t") for f in absolute_segments],
                   ignore_index=True)[abscolumns]
    assert(len(absdf) == len(absegdf))
    for c in abscolumns:
        absdf[c] = absegdf[c]
    absdf["clonality"] = absdf.apply(lambda x: "clonal" if x.Pr_somatic_clonal >= .5 or x.ccf_CI95_low >= .9 else "subclonal", axis=1)

    return absdf


def add_loh(df):
    rv = df.copy()
    if len(df) > 0:
        if "facetsLCN_EM" in df.columns:
            rv["LOH"] = df.facetsLCN_EM.apply(lambda x: "true" if x == 0 or x =="0" else ".")
        else:
            rv["LOH"] = pd.Series(["N/A"] * len(df))
    else:
        rv["LOH"] = pd.Series()

    return rv


def add_maf(df):
    rv = df.copy()
    def f(x):
        if "," in x["TUMOR.AD"]:
            return float(x["TUMOR.AD"].split(",")[1]) / sum(map(float, x["TUMOR.AD"].split(',')))
        else:
            return float('nan')
    def g(x):
        if "," in x["NORMAL.AD"]:
            return float(x["NORMAL.AD"].split(",")[1]) / sum(map(float, x["NORMAL.AD"].split(',')))
        else:
            return float('nan')
    if len(df) > 0:
        rv["TUMOR_MAF"] = df.apply(f, axis=1)
        rv["NORMAL_MAF"] = df.apply(g, axis=1)
    else:
        rv["TUMOR_MAF"] = pd.Series()
        rv["NORMAL_MAF"] = pd.Series()
    return rv



def add_cancer_gene(df):
    rv = df.copy()
    rv["cancer_gene"] = df.apply(lambda x: "true" if x["cancer_gene_census"] == "true" or
                                 x["kandoth"] == "true" or
                                 x["lawrence"] == "true" else ".",
                                 axis=1)
    return rv


def add_non_existent_columns(df, columns, fill_value):
    rv = df.copy()
    for c in columns:
        if c not in df.columns:
            rv[c] = len(df) * [fill_value]
    return rv


def add_columns_write_excel(df, writer, sheetname, absdf=None, write_columns=None, output_tsv_dir=None, annotdf=None):
    if all([c in df.columns for c in "TUMOR.AD NORMAL.AD".split()]):
        df = add_maf(df)
    if len(df > 0):
        if all([c in df.columns for c in "cancer_gene_census kandoth lawrence".split()]):
            df = add_cancer_gene(df)
        df = add_loh(df)
        if write_columns:
            df = df[[c for c in write_columns if c in df.columns]]
        df = df.set_index("TUMOR_SAMPLE NORMAL_SAMPLE CHROM POS REF ALT".split())
        if absdf is not None:
            df = df.join(absdf.set_index("TUMOR_SAMPLE NORMAL_SAMPLE CHROM POS REF ALT".split()), how='left')
        if annotdf is not None:
            df = df.join(annotdf.set_index("TUMOR_SAMPLE NORMAL_SAMPLE CHROM POS REF ALT".split()), how='left')
        df.reset_index().to_excel(writer, sheetname, index=False)
        if output_tsv_dir:
            df.to_csv(output_tsv_dir + "/" + sheetname.lower() + ".tsv", sep="\t", index=True)


def write_mutation_summary(mutect_snps_high_moderate, mutect_snps_low_modifier,
                           mutect_snps_synonymous, mutect_snps_nonsynonymous,
                           mutect_indels_high_moderate,
                           mutect_indels_low_modifier,
                           mutect_indels_synonymous, mutect_indels_nonsynonymous,
                           hotspot,
                           excel_file,
                           absolute_somatic_txts,
                           absolute_segments,
                           output_tsv_dir,
                           annotation_tsv,
                           max_exac_af):
    # create output tsv dir if required
    if output_tsv_dir:
        mkdir_p(output_tsv_dir)

    # create absolute df with cancer cell fractions
    if absolute_somatic_txts and absolute_segments:
        absdf = create_absolute_df(absolute_somatic_txts, absolute_segments)
    else:
        absdf = None
    # create annotation df
    if annotation_tsv:
        annotdf = pd.read_csv(annotation_tsv, sep="\t")
    else:
        annotdf = None
    summary_columns = "CHROM,POS,TUMOR_SAMPLE,NORMAL_SAMPLE,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].HGVS_C,ANN[*].EFFECT,TUMOR_MAF,NORMAL_MAF,TUMOR.DP,NORMAL.DP,ExAC_AF,dbNSFP_MutationTaster_pred,fathmm_pred,dbNSFP_PROVEAN_pred,LOH,pathogenicity,HOTSPOT".split(",")
    # find chasm score columns, they are prefixed with chosen classifier
    chasm_score_columns = [c for c in pd.read_csv(mutect_snps_high_moderate, sep="\t").columns if "chasm_score" in c]
    # add gene annotations and chasm score columns
    summary_columns += chasm_score_columns + "cancer_gene_census,kandoth,lawrence,hap_insuf,REF,ALT,ANN[*].IMPACT".split(",")

    writer = pd.ExcelWriter(excel_file)

    def read_tsv(tsv):
        return pd.read_csv(tsv, sep="\t", dtype={"CHROM": str})

    # add summaries
    required_columns = summary_columns + "NORMAL.AD TUMOR.AD facetsLCN_EM".split()
    mutsdf = pd.concat([add_non_existent_columns(filter_annotations_with_impact(read_tsv(mutect_snps_high_moderate), "HIGH|MODERATE"), required_columns, ".")[required_columns],
                        add_non_existent_columns(filter_annotations_with_impact(read_tsv(mutect_snps_low_modifier), "LOW", effect=".*synonymous_variant.*"), required_columns, ".")[required_columns],
                        add_non_existent_columns(filter_annotations_with_impact(read_tsv(mutect_indels_high_moderate), "HIGH|MODERATE"), required_columns, ".")[required_columns]],
                        ignore_index=True).sort_values("TUMOR_SAMPLE CHROM POS".split())
    exac_af_sel = mutsdf["ExAC_AF"].apply(lambda x: x == "." or float(x) < max_exac_af)
    add_columns_write_excel(mutsdf[exac_af_sel], writer, "MUTATION_SUMMARY", absdf,
        write_columns=summary_columns, output_tsv_dir=output_tsv_dir,
        annotdf=annotdf)
    add_columns_write_excel(filter_annotations_with_impact(read_tsv(mutect_snps_high_moderate), "HIGH|MODERATE"),
        writer, "SNV_HIGH_MODERATE_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_snps_low_modifier), writer, "SNV_LOW_MODIFIER_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_snps_synonymous), writer, "SNV_SYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(filter_annotations_with_impact(read_tsv(mutect_snps_nonsynonymous), "HIGH|MODERATE"),
        writer, "SNV_NONSYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(filter_annotations_with_impact(read_tsv(mutect_indels_high_moderate), "HIGH|MODERATE"),
        writer, "INDEL_HIGH_MODERATE_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_indels_low_modifier), writer, "INDEL_LOW_MODIFIER_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_indels_synonymous), writer, "INDEL_SYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(filter_annotations_with_impact(read_tsv(mutect_indels_nonsynonymous), "HIGH|MODERATE"),
        writer, "INDEL_NONSYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns, output_tsv_dir=output_tsv_dir, annotdf=annotdf)

    # add raw files both as excel and tsv
    add_columns_write_excel(read_tsv(hotspot), writer, "hotspot", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_snps_high_moderate), writer, "mutect_snps_high_moderate", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_snps_low_modifier), writer, "mutect_snps_low_modifier", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_snps_synonymous), writer, "mutect_snps_synonymous", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_snps_nonsynonymous), writer, "mutect_snps_nonsynonymous", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_indels_high_moderate), writer, "mutect_indels_high_moderate", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_indels_low_modifier), writer, "mutect_indels_low_modifier", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_indels_synonymous), writer, "mutect_indels_synonymous", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)
    add_columns_write_excel(read_tsv(mutect_indels_nonsynonymous), writer, "mutect_indels_nonsynonymous", absdf, output_tsv_dir=output_tsv_dir, annotdf=annotdf)

    writer.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mutect_snps_high_moderate", type=str, help="TSV")
    parser.add_argument("mutect_snps_low_modifier", type=str, help="TSV")
    parser.add_argument("mutect_snps_synonymous", type=str, help="TSV")
    parser.add_argument("mutect_snps_nonsynonymous", type=str, help="TSV")
    parser.add_argument("mutect_indels_high_moderate", type=str, help="TSV")
    parser.add_argument("mutect_indels_low_modifier", type=str, help="TSV")
    parser.add_argument("mutect_indels_synonymous", type=str, help="TSV")
    parser.add_argument("mutect_indels_nonsynonymous", type=str, help="TSV")
    parser.add_argument("hotspot", type=str, help="TSV")
    parser.add_argument("excel_file", type=str, help="mutation summary excel")
    parser.add_argument("--absolute_somatic_txts", default=None, type=str, help="TSV comma separated list of somatic files of absolute input")
    parser.add_argument("--absolute_segments", default=None, type=str, help="TSV comma separated list of absolute mutations output")
    parser.add_argument("--output_tsv_dir", default=None, type=str, help="Output raw sheets as tsv in given directory")
    parser.add_argument("--annotation_tsv", default=None, type=str, help="File with TUMOR_SAMPLE NORMAL_SAMPLE CHROM POS REF ALT, plus other columns of choice for annotation")
    parser.add_argument("--max_exac_af", default=1, type=float, help="Set threshold for ExAC_AF column. Only applied to MUTATION_SUMMARY column.")
    args = parser.parse_args()
    if args.absolute_somatic_txts and args.absolute_segments:
        absolute_somatic_txts = args.absolute_somatic_txts.split(",")
        absolute_segments = args.absolute_segments.split(",")
        if not (len(absolute_somatic_txts) == len(absolute_segments)):
            raise(Exception("Unequal length of absolute files"))
    else:
        absolute_somatic_txts = None
        absolute_segments = None
    write_mutation_summary(args.mutect_snps_high_moderate,
                           args.mutect_snps_low_modifier,
                           args.mutect_snps_synonymous,
                           args.mutect_snps_nonsynonymous,
                           args.mutect_indels_high_moderate,
                           args.mutect_indels_low_modifier,
                           args.mutect_indels_synonymous,
                           args.mutect_indels_nonsynonymous,
                           args.hotspot,
                           args.excel_file,
                           absolute_somatic_txts,
                           absolute_segments,
                           args.output_tsv_dir,
                           args.annotation_tsv,
                           args.max_exac_af)

if __name__ == "__main__":
    main()
