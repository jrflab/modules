#!/usr/bin/env python

""" classify provean of vcf records
"""

import vcf
import argparse
import sys
import remote_provean_query as rpq

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='provean_vcf.py',
                                     description='Add provean to vcf file')
    parser.add_argument('vcf_infile')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))
    vcf_reader.infos['provean_protein_id'] = vcf.parser._Info(id='provean_protein_id', num=1, type='String',
                                                              desc="provean protein id",
                                                              source=None, version=None)
    vcf_reader.infos['provean_pred'] = vcf.parser._Info(id='provean_pred', num=1, type='String',
                                                        desc="Mutation taster prediction using webquery if indel",
                                                        source=None, version=None)
    vcf_reader.infos['provean_score'] = vcf.parser._Info(id='provean_score', num=1, type='Float',
                                                         desc="Mutation taster score using webquery if indel",
                                                         source=None, version=None)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    records = [record for record in vcf_reader]
    query_records = []

    for record in records:
        if record.is_indel:
            query_records.append(record)

    if len(query_records) > 0:
        query = rpq.RemoteProveanQuery(query_records)
        query.run_query()

    for record in records:
        vcf_writer.write_record(record)
    vcf_writer.close()
