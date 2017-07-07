#!/usr/bin/env python

import glob2
import yaml
import argparse
import re

parser = argparse.ArgumentParser(prog='create_iontorrent_sample_yaml.py',
                                 description='Create samples.merge.yaml and samples.yaml(best guess) from a rawdata \
                                 dir')
parser.add_argument('sample_map')
parser.add_argument('bam_dir', nargs='+')
parser.add_argument('--sample_merge_file', help='sample fastq file yaml output',
                    type=argparse.FileType('w', 0), nargs='?',
                    default='sample.merge.yaml')
parser.add_argument('--sample_file', help='sample yaml output file', type=argparse.FileType('w', 0), nargs='?',
                    default='samples.yaml')
args = parser.parse_args()

sample_map = {}
with open(args.sample_map, 'r') as f:
    for line in f:
        sline = line.split()
        s = sline[0]
        if s not in sample_map:
            sample_map[s] = []
        sample_map[sline[0]].append(sline[1:])
f.close()

bam_files = []
for d in args.bam_dir:
    bam_files.extend(glob2.glob(d + '/**/*.bam'))

sample_merge = {}
for sample, keys in list(sample_map.items()):
    if sample not in sample_merge:
        sample_merge[sample] = []
    for key in keys:
        for bam_file in bam_files:
            if all([k in bam_file for k in key]):
                sample_merge[sample].append(bam_file)

yaml.dump(sample_merge, args.sample_merge_file)

# output best guess for samples.yaml
normals = set()
tumors = set()
for s in list(sample_merge.keys()):
    if s.upper().endswith('N'):
        normals.add(s)
    else:
        tumors.add(s)

samples = []
unmatched_tumors = tumors
for norm in normals:
    name = re.sub(r'[nN]$', '', norm)
    tums = [x for x in tumors if x.startswith(name)]
    if tums is not None:
        for x in tums:
            unmatched_tumors.discard(x)
        samples.append({'name': name, 'normal': norm, 'tumor': tums})
    else:
        samples.append({'name': name, 'normal': norm})

for unmatched_tum in unmatched_tumors:
    name = re.sub(r'[tT]$', '', unmatched_tum)
    samples.append({'name': name, 'tumor': [unmatched_tum]})

yaml.dump(samples, args.sample_file)
