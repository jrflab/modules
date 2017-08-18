#!/usr/bin/env python
"""
Annotate genotype to sufam vcf file
"""

import argparse
import vcf
import collections
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sufam_vcf_file', help='multi-sample sufam file')
    parser.add_argument('vcf_files', nargs='+', help='sample pair vcf files')
    args = parser.parse_args()

    sample_variants = collections.defaultdict(set)
    for f in args.vcf_files:
        vcf_reader = vcf.Reader(open(f, 'r'))
        for record in vcf_reader:
            recid = "{}:{}:{}/{}".format(record.CHROM, record.POS, record.REF, record.ALT)
            s = record.samples[0].sample
            sample_variants[s].add(recid)

    sufam_vcf_reader = vcf.Reader(open(args.sufam_vcf_file, 'r'))
    sufam_vcf_reader.infos['samples_called_in'] = vcf.parser._Info(id='samples_called_in', num='.',
                                                                   type='String',
                                                                   desc='samples called in',
                                                                   source=None,
                                                                   version=None)
    vcf_writer = vcf.Writer(sys.stdout, sufam_vcf_reader)
    for record in sufam_vcf_reader:
        recid = "{}:{}:{}/{}".format(record.CHROM, record.POS, record.REF, record.ALT)
        for s, v in list(sample_variants.items()):
            if recid in v:
                if 'samples_called_in' not in record.INFO:
                    record.INFO['samples_called_in'] = []
                record.INFO['samples_called_in'].append(s)
        vcf_writer.write_record(record)
    vcf_writer.close()
