#!/usr/bin/env python

import glob2
import yaml
import argparse
import re

parser = argparse.ArgumentParser(prog='create_sample_yaml.py',
                                 description='Create samples.fastq.yaml and samples.yaml(best guess) from a rawdata \
                                 dir')
parser.add_argument('fastq_dir', nargs = '+')
parser.add_argument('--fastq_suffix', default='.fastq.gz')
parser.add_argument('--sample_fastq_file', help='sample fastq file yaml output',
                    type=argparse.FileType('w'), nargs='?',
                    default='sample.fastq.yaml')
parser.add_argument('--sample_file', help='sample yaml output file', type=argparse.FileType('w'), nargs='?',
                    default='samples.yaml')
args = parser.parse_args()

paired_fastqs = []
for fastq_dir in args.fastq_dir:
    fastqFiles = glob2.glob(fastq_dir + '/**/*' + args.fastq_suffix)
    r1fastqs = [x for x in fastqFiles if re.search(r'_S\d+_R1_', x)]
    for r1fastq in r1fastqs:
        r2fastq = re.sub('_(S\d+)_R1_', '_\g<1>_R2_', r1fastq)
        assert(any([r2fastq == x for x in fastqFiles]))  # r2 fastq doesn't exist
        paired_fastqs.append([r1fastq, r2fastq])

sample_fastqs = {}
for pair in paired_fastqs:
    sample = re.search(r'Sample_([^/]+)', pair[0]).group(1)
    sample = re.sub('_IGO.*', '', sample)
    sample = re.sub('_', '-', sample)
    if sample not in sample_fastqs:
        sample_fastqs[sample] = []
    sample_fastqs[sample].append(pair)

yaml.dump(sample_fastqs, args.sample_fastq_file)

# output best guess for samples.yaml
normals = set()
tumors = set()
for s in list(sample_fastqs.keys()):
    if s.endswith('N'):
        normals.add(s)
    else:
        tumors.add(s)

samples = []
unmatched_tumors = tumors
for norm in normals:
    name = re.sub(r'N$', '', norm)
    tums = [x for x in tumors if x.startswith(name)]
    if tums is not None:
        for x in tums:
            unmatched_tumors.discard(x)
        samples.append({'name': name, 'normal': norm, 'tumor': tums})
    else:
        samples.append({'name': name, 'normal': norm})

for unmatched_tum in unmatched_tumors:
    name = re.sub(r'T$', '', unmatched_tum)
    samples.append({'name': name, 'tumor': [unmatched_tum]})

yaml.dump(samples, args.sample_file)
