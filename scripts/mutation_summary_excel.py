#!/usr/bin/env python
"""
Create mutation summary excel from given tsvs
"""
import argparse
import pandas as pd
import os


def create_absolute_df(absolute_somatic_txts, absolute_segments):
    absdf = pd.concat([pd.read_csv(asegs, sep="\t", dtype={"Chromosome":str}) for asegs
               in absolute_somatic_txts], ignore_index=True)
    absdf["Chromosome"] = absdf.Chromosome.replace({"23":"X","24":"Y"})
    absdf["TUMOR_SAMPLE"] = absdf.Sample.apply(lambda x: x.split("_")[0])
    absdf["NORMAL_SAMPLE"] = absdf.Sample.apply(lambda x: x.split("_")[1])
    absdf.drop([c for c in absdf.columns if c not in "TUMOR_SAMPLE NORMAL_SAMPLE Ref Alt Chromosome Position".split()],
                         axis=1,
                         inplace=True)
    absdf.rename(columns={"Ref":"REF","Alt":"ALT","Chromosome":"CHROM","Position":"POS"}, inplace=True)
    ccfrac = pd.concat([pd.read_csv(f, dtype={"Chromosome":str}, sep="\t") for f in absolute_segments],
                   ignore_index=True)["cancer_cell_frac"]
    assert(len(absdf) == len(ccfrac))
    absdf["cancer_cell_frac"] = ccfrac

    return absdf

def add_maf(df):
    if len(df) > 0:
        df["TUMOR_MAF"] = df.apply(lambda x: float(x["TUMOR.AD"].split(",")[1]) / x["TUMOR.DP"], axis=1)
        df["NORMAL_MAF"] = df.apply(lambda x: float(x["NORMAL.AD"].split(",")[1]) / x["NORMAL.DP"], axis=1)
    else:
        df["TUMOR_MAF"] = pd.Series()
        df["NORMAL_MAF"] = pd.Series()
    return df


def add_maf_join_abs_write_excel(tsv, writer, sheetname, absdf=None, write_columns=None):
    df = add_maf(pd.read_csv(tsv, sep="\t", dtype={"CHROM":str}))
    if len(df > 0):
        if write_columns:
            df = df[[c for c in write_columns if c in df.columns]]
        if absdf is not None:
            df = df.set_index("TUMOR_SAMPLE NORMAL_SAMPLE CHROM POS REF ALT".split())\
                .join(absdf.set_index("TUMOR_SAMPLE NORMAL_SAMPLE CHROM POS REF ALT".split()), how='left')
        df.to_excel(writer, sheetname, index=(absdf is not None))


def write_mutation_summary(mutect_high_moderate, mutect_low_modifier,
                           mutect_synonymous, mutect_nonsynonymous,
                           strelka_varscan_high_moderate,
                           strelka_varscan_low_modifier,
                           strelka_varscan_synonymous, strelka_varscan_nonsynonymous,
                           excel_file,
                           absolute_somatic_txts=None,
                           absolute_segments=None):
    # create absolute df with cancer cell fractions
    if absolute_somatic_txts and absolute_segments:
        absdf = create_absolute_df(absolute_somatic_txts, absolute_segments)
    else:
        absdf = None
    summary_columns = "CHROM,POS,TUMOR_SAMPLE,NORMAL_SAMPLE,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].EFFECT,TUMOR_MAF,NORMAL_MAF,TUMOR.DP,NORMAL.DP,dbNSFP_MutationTaster_pred,fathmm_pred".split(",")
    # find chasm score columns, they are prefixed with chosen classifier
    chasm_score_columns = [c for c in pd.read_csv(mutect_high_moderate, sep="\t").columns if "chasm_score" in c]
    # add gene annotations and chasm score columns
    summary_columns += chasm_score_columns + "cancer_gene_census,kandoth,lawrence,hap_insuf,REF,ALT,ANN[*].IMPACT".split(",")

    writer = pd.ExcelWriter(excel_file)

    # add summaries
    add_maf_join_abs_write_excel(mutect_high_moderate, writer, "SNV_HIGH_MODERATE_SUMMARY", absdf, write_columns=summary_columns)
    add_maf_join_abs_write_excel(mutect_low_modifier, writer, "SNV_LOW_MODIFIER_SUMMARY", absdf, write_columns=summary_columns)
    add_maf_join_abs_write_excel(mutect_synonymous, writer, "SNV_SYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns)
    add_maf_join_abs_write_excel(mutect_nonsynonymous, writer, "SNV_NONSYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns)
    add_maf_join_abs_write_excel(strelka_varscan_high_moderate, writer, "INDEL_HIGH_MODERATE_SUMMARY", absdf, write_columns=summary_columns)
    add_maf_join_abs_write_excel(strelka_varscan_low_modifier, writer, "INDEL_LOW_MODIFIER_SUMMARY", absdf, write_columns=summary_columns)
    add_maf_join_abs_write_excel(strelka_varscan_synonymous, writer, "INDEL_SYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns)
    add_maf_join_abs_write_excel(strelka_varscan_nonsynonymous, writer, "INDEL_NONSYNONYMOUS_SUMMARY", absdf, write_columns=summary_columns)

    # add raw files
    add_maf_join_abs_write_excel(mutect_high_moderate, writer, "mutect_high_moderate", absdf)
    add_maf_join_abs_write_excel(mutect_low_modifier, writer, "mutect_low_modifier", absdf)
    add_maf_join_abs_write_excel(strelka_varscan_high_moderate, writer, "strelka_varscan_high_moderate", absdf)
    add_maf_join_abs_write_excel(strelka_varscan_low_modifier, writer, "strelka_varscan_low_modifier", absdf)
    add_maf_join_abs_write_excel(strelka_varscan_synonymous, writer, "strelka_varscan_synonymous", absdf)
    add_maf_join_abs_write_excel(strelka_varscan_nonsynonymous, writer, "strelka_varscan_nonsynonymous", absdf)

    writer.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mutect_high_moderate", type=str, help="TSV")
    parser.add_argument("mutect_low_modifier", type=str, help="TSV")
    parser.add_argument("mutect_synonymous", type=str, help="TSV")
    parser.add_argument("mutect_nonsynonymous", type=str, help="TSV")
    parser.add_argument("strelka_varscan_high_moderate", type=str, help="TSV")
    parser.add_argument("strelka_varscan_low_modifier", type=str, help="TSV")
    parser.add_argument("strelka_varscan_synonymous", type=str, help="TSV")
    parser.add_argument("strelka_varscan_nonsynonymous", type=str, help="TSV")
    parser.add_argument("excel_file", type=str, help="mutation summary excel")
    parser.add_argument("--absolute_somatic_txts", default=None, type=str, help="TSV comma separated list of somatic files of absolute input")
    parser.add_argument("--absolute_segments", default=None, type=str, help="TSV comma separated list of absolute mutations output")
    args = parser.parse_args()
    if args.absolute_somatic_txts and args.absolute_segments:
        absolute_somatic_txts = args.absolute_somatic_txts.split(",")
        absolute_segments = args.absolute_segments.split(",")
        if not (len(absolute_somatic_txts) == len(absolute_segments)):
            raise(Exception("Unequal length of absolute files"))
    else:
        absolute_somatic_txts = None
        absolute_segments = None
    write_mutation_summary(args.mutect_high_moderate,
                           args.mutect_low_modifier,
                           args.mutect_synonymous,
                           args.mutect_nonsynonymous,
                           args.strelka_varscan_high_moderate,
                           args.strelka_varscan_low_modifier,
                           args.strelka_varscan_synonymous,
                           args.strelka_varscan_nonsynonymous,
                           args.excel_file,
                           absolute_somatic_txts,
                           absolute_segments)

if __name__ == "__main__":
    main()
