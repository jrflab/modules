#!/usr/bin/env python

import argparse
import math
import vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='split_vcf.py',
                                     description='split a vcf file into chunks')
    parser.add_argument('--num_chunks', '-c', type=int, default=100, help='number of chunks')
    parser.add_argument('--out_prefix', '-o', required=True, help='output prefix')
    parser.add_argument('vcf_file', help='bed file to split')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_file, 'r'))
    records = [record for record in vcf_reader]

    n = int(math.ceil(len(records) / args.num_chunks))
    x = 1
    for i in range(0, len(records), n):
        vcf_writer = vcf.Writer(open(args.out_prefix + '{0:03d}.vcf'.format(x), 'w'), vcf_reader)
        if x == args.num_chunks:
            # last chunk, write everything
            for record in records[i:]:
                vcf_writer.write_record(record)
            vcf_writer.close()
            break
        else:
            for record in records[i:i + n]:
                vcf_writer.write_record(record)
            vcf_writer.close()
        x = x + 1
