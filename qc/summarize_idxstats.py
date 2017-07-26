#!/usr/bin/env python

""" collate the samtools idxstats outputs"""

import argparse
import pycurl
import re
import os
import sys
import pandas as pd
from functools import reduce


def slack_warning(slack_group, slack_token, slack_channel, message):
    slackbot_url = "https://{}.slack.com/services/hooks/slackbot" \
        "?token={}&channel=%23{}".format(slack_group, slack_token, slack_channel)
    c = pycurl.Curl()
    c.setopt(c.URL, slackbot_url)
    c.setopt(c.POSTFIELDS, message)
    c.setopt(pycurl.WRITEFUNCTION, lambda x: None)
    c.perform()
    c.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='summarize_idxstats.py',
                                     description='summarize samtools idxstats; alert if missing coverage for contig')
    parser.add_argument('idxstats_file', nargs='+')
    parser.add_argument('--targets_file', nargs='?', help='targets bed file')
    parser.add_argument('--excel_file', nargs='?', help='excel file to write')
    parser.add_argument('--project_name', nargs='?', help='name of the project (for notifications)')
    parser.add_argument('--slack_group', nargs='?', default='jrflab', help='slack group for notifications')
    parser.add_argument('--slack_token', nargs='?', default='x2UXxJD2PlB1VBurY5MXIlW8', help='slack token')
    parser.add_argument('--slack_channel', nargs='?', default='pipeline_qc', help='slack channel for notifications')
    parser.add_argument('--chromosomes', nargs='*', help='chromosomes')
    args = parser.parse_args()

    idxs = []
    for idxstats_file in args.idxstats_file:
        s = os.path.splitext(os.path.basename(idxstats_file))[0]
        df = pd.read_table(idxstats_file, sep='\t', header=None)
        df.columns = ['chrom', 'length', s + '.num_mapped', s + '.num_unmapped']
        idxs.append(df)

    idxstats = reduce(lambda x, y: x.merge(y, how='inner', on=['chrom', 'length']), idxs)

    if args.targets_file is not None:
        targets = pd.read_table(args.targets_file, sep='\t', header=None, dtype={0: str})
        contigs = targets[0].unique().tolist()
    elif args.chromosomes is not None:
        contigs = args.chromosomes
    else:
        contigs = idxstats.chrom.unique().tolist()

    contigs = [x for x in contigs if 'Y' not in x]

    num_mapped = idxstats.loc[idxstats['chrom'].isin(contigs)].select(lambda x: x.endswith('num_mapped'), axis=1)
    missing_contig_cov_samples = []
    for (col, v) in list(num_mapped.items()):
        s = re.sub(r'\..*', r'', col)
        if any(v == 0):
            missing_contig_cov_samples.append(s)

    if len(missing_contig_cov_samples) > 0:
        samples = ' '.join(missing_contig_cov_samples)
        warning_message = ":suspect: missing contigs : {} : {}".format(samples, os.getcwd())
        if args.project_name is not None:
            warning_message = args.project_name + ": " + warning_message
        slack_warning(args.slack_group, args.slack_token, args.slack_channel, warning_message)

    idxstats.to_csv(sys.stdout, sep='\t', index=False)

    if args.excel_file is not None:
        writer = pd.ExcelWriter(args.excel_file)
        idxstats.to_excel(writer, 'idxstats', index=False)
        writer.close()
