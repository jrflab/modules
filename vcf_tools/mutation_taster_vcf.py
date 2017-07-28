#!/usr/bin/env python

import argparse
import vcf
import pandas as pd
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='mutation_taster_vcf.py',
                                     description='Add mutation taster results to vcf file')
    parser.add_argument('vcf_infile')

    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    vcf_reader.infos['MutationTaster_pred'] = vcf.parser._Info(id='MutationTaster_pred', num=1, type='String',
                                                            desc="Mutation taster prediction using webquery if indel or"
                                                            " dbNSFP if it exists", source=None, version=None)
    vcf_reader.infos['MutationTaster_score'] = vcf.parser._Info(id='MutationTaster_score', num=1, type='Float',
                                                                desc="Mutation taster score using webquery if indel or"
                                                                " dbNSFP if it exists", source=None, version=None)

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    prediction_priority = ['disease_causing_automatic',
                        'disease_causing',
                        'polymorphism_automatic',
                        'polymorphism']

    for record in vcf_reader:
        pred = None
        score = None
        if ('dbNSFP_MutationTaster_pred' in record.INFO and record.INFO['dbNSFP_MutationTaster_pred'][0] is not None):
            pred = record.INFO['dbNSFP_MutationTaster_pred'][0]
            score = record.INFO['dbNSFP_MutationTaster_score'][0]
        elif record.is_indel:
            url = "http://www.mutationtaster.org/cgi-bin/MutationTaster/MT_ChrPos.cgi" \
                "?chromosome={}&position={}&ref={}&alt={}".format(record.CHROM, record.POS, record.REF, record.ALT[0])
            sys.stderr.write("Querying: {}\n".format(url))
            dfs = pd.read_html(url)
            if (len(dfs) < 2):
                sys.stderr.write("query failed, skipping\n")
            else:
                # summary is the second dataframe
                summary = dfs[1][1:]
                summary.columns = dfs[1].iloc[0]
                if ('prediction' in summary.columns):
                    pred_score = {}
                    for i, row in summary.iterrows():
                        if (row['prediction'] not in pred_score):
                            pred_score[row['prediction']] = []
                        pred_score[row['prediction']].append(float(row['probability']))
                    for p in prediction_priority:
                        if p in pred_score:
                            pred = p
                            score = max(pred_score[p])
                            break
        if pred is not None and score is not None:
            record.INFO['MutationTaster_pred'] = pred
            record.INFO['MutationTaster_score'] = score
        vcf_writer.write_record(record)
    vcf_writer.close()
