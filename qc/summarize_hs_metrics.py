#!/usr/bin/env python

"""summarize_hs_metrics.py: summarize the hs metrics file"""

import sys
import os
import collections
import argparse
import pycurl
import pandas as pd

column_map = collections.OrderedDict([
    ('Sample ID', 'SAMPLE'),
    ('Target Territory', 'TARGET_TERRITORY'),
    ('Total Reads', 'TOTAL_READS'),
    ('Percent Selected Bases', 'PCT_SELECTED_BASES'),
    ('Mean Target Coverage', 'MEAN_TARGET_COVERAGE'),
    ('Percent Target Bases 2X', 'PCT_TARGET_BASES_2X'),
    ('Percent Target Bases 10X', 'PCT_TARGET_BASES_10X'),
    ('Percent Target Bases 20X', 'PCT_TARGET_BASES_20X'),
    ('Percent Target Bases 30X', 'PCT_TARGET_BASES_30X'),
    ('Percent Target Bases 40X', 'PCT_TARGET_BASES_40X'),
    ('Percent Target Bases 50X', 'PCT_TARGET_BASES_50X'),
    ('Percent Target Bases 100X', 'PCT_TARGET_BASES_100X')
])

def slack_warning(slack_group, slack_token, slack_channel, message):
    slackbot_url = "https://{}.slack.com/services/hooks/slackbot" \
        "?token={}&channel=%23{}".format(slack_group, slack_token,slack_channel)
    c = pycurl.Curl()
    c.setopt(c.URL, slackbot_url)
    c.setopt(c.POSTFIELDS, message)
    c.setopt(pycurl.WRITEFUNCTION, lambda x: None)
    c.perform()
    c.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='summarize_hs_metrics.py',
                                     description='summarize HS metrics; alert users if is cov below a threshold')
    parser.add_argument('hs_metrics_file', nargs='+')
    parser.add_argument('--depth_threshold', nargs='?', default=100, help='depth threshold for sending alerts')
    parser.add_argument('--project_name', nargs='?', help='name of the project (for notifications)')
    parser.add_argument('--slack_group', nargs='?', default='jrflab', help='slack group for notifications')
    parser.add_argument('--slack_token', nargs='?', default='x2UXxJD2PlB1VBurY5MXIlW8', help='slack token')
    parser.add_argument('--slack_channel', nargs='?', default='pipeline_qc', help='slack channel for notifications')
    parser.add_argument('--excel_file', nargs='?', help='excel file to write')
    args = parser.parse_args()

    metrics = pd.DataFrame()
    for metrics_file in args.hs_metrics_file:
        sample_metrics = pd.read_table(metrics_file, comment='#')
        sample_name = os.path.basename(metrics_file).split('.')[0] # remove suffix
        sample_metrics = sample_metrics[column_map.values()].copy()
        sample_metrics.loc[0, 'SAMPLE'] = sample_name
        metrics = metrics.append(sample_metrics)

    below_threshold = metrics.query('MEAN_TARGET_COVERAGE < @args.depth_threshold')
    if len(below_threshold) > 0:
        samples = ' '.join(below_threshold['SAMPLE'])
        warning_message = ":suspect: <{}X cov: {} : {}".format(args.depth_threshold, samples, os.getcwd())
        if args.project_name is not None:
            warning_message = args.project_name + ": " + warning_message
        slack_warning(args.slack_group, args.slack_token, args.slack_channel, warning_message)

    metrics.columns = column_map.keys()

    metrics.to_csv(sys.stdout, sep='\t', index=False)

    if args.excel_file is not None:
        writer = pd.ExcelWriter(args.excel_file)
        metrics.to_excel(writer, 'Sequencing Statistics', index=False)
        writer.close()
