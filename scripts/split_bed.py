#!/usr/bin/env python

import argparse
import math

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='split_bed.py',
                                     description='split a bed file into chunks')
    parser.add_argument('--num_chunks', '-c', type=int, default=100, help='number of chunks')
    parser.add_argument('--out_prefix', '-o', required=True, help='output prefix')
    parser.add_argument('bed_file', help='bed file to split')
    args = parser.parse_args()

    bed = [line for line in open(args.bed_file, 'r')]

    n = int(math.ceil(len(bed) / args.num_chunks))
    x = 1
    for i in range(0, len(bed), n):
        f = open(args.out_prefix + '{0:03d}.bed'.format(x), 'w')
        if x == args.num_chunks:
            # last chunk, write everything
            for line in bed[i:]:
                f.write(line)
            f.close()
            break
        else:
            for line in bed[i:i + n]:
                f.write(line)
            f.close()
        x = x + 1
