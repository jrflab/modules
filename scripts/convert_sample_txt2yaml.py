#!/usr/bin/env python

import yaml
import argparse

parser = argparse.ArgumentParser(prog='convert_sample_txt2yaml.py',
                                 description='Convert samples.txt/sample_sets.txt to yaml')
parser.add_argument('sample_txt_file')
parser.add_argument('--out_file', help='sample yaml output file', type=argparse.FileType('w', 0), nargs='?',
                    default='samples.yaml')
args = parser.parse_args()

samples = []
with open(args.sample_txt_file, 'r') as f:
    for line in f:
        split_sp = line.rstrip().split()
        if (len(split_sp) == 1):
            samples.append({'tumor': split_sp})
        else:
            samples.append({'normal': split_sp[-1], 'tumor': split_sp[:-1]})

yaml.dump(samples, args.out_file)
