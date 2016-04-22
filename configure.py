#!/usr/bin/env python
from __future__ import print_function
import yaml
import argparse

parser = argparse.ArgumentParser(prog='configure',
                                 description='Convert project YAML file to Make')
parser.add_argument('--project_config_file', type=file, help='project yaml config file',
                    default='project_config.yaml')
parser.add_argument('--config_file', type=file, help='modules yaml config file',
                    default='modules/config.yaml')
parser.add_argument('--samples_file', type=file, help='yaml samples file', default='samples.yaml')
parser.add_argument('--sample_fastq_file', help='yaml sample fastq file mappings', type=file,
                    default='sample.fastq.yaml')
parser.add_argument('--sample_merge_file', help='yaml sample merge mappings', type=file)
parser.add_argument('--outFile', help='project make include file', type=argparse.FileType('w', 0), nargs='?',
                    default='project_config.inc')
args = parser.parse_args()


def lowerBool(x):
    if isinstance(x, bool):
        return str(x).lower()
    else:
        return x

of = args.outFile

config = yaml.load(args.project_config_file)
for k, v in config.iteritems():
    print("{} = {}".format(k.upper(), lowerBool(v)), file=of)

samples = yaml.load(args.samples_file)

tumors = set()
normals = set()
normal_tumors = {}
tumor_normal = {}
for s in samples:
    if 'normal' in s:
        normals.add(s['normal'])
    if 'tumor' in s:
        for t in s['tumor']:
            tumors.add(t)
    if 'normal' in s and 'tumor' in s:
        normal_tumors[s['normal']] = s['tumor']
        for t in s['tumor']:
            tumor_normal[t] = s['normal']

print("", file=of)
print("SAMPLES = {} {}".format(" ".join(tumors), " ".join(normals)), file=of)
print("TUMOR_SAMPLES = {}".format(" ".join(tumors)), file=of)
print("NORMAL_SAMPLES = {}".format(" ".join(normals)), file=of)
print("SAMPLE_PAIRS = ", end="", file=of)
for k, v in tumor_normal.iteritems():
    print("{}_{}".format(k, v), end=" ", file=of)
print("", file=of)
print("SAMPLE_SETS = ", end="", file=of)
for k, v in normal_tumors.iteritems():
    print("{}_{}".format("_".join(v), k), end=" ", file=of)
print("", file=of)

print("\n# tumor.normal = tumors", file=of)
for k, v in normal_tumors.iteritems():
    print("tumor.{} = {}".format(k, " ".join(v)), file=of)
print("\n# normal.tumor = normal", file=of)
for k, v in tumor_normal.iteritems():
    print("normal.{} = {}".format(k, v), file=of)
print("\n# tumor.pair = tumor", file=of)
for k, v in tumor_normal.iteritems():
    print("tumor.{}_{} = {}".format(k, v, k), file=of)
print("\n# normal.pair = normal", file=of)
for k, v in tumor_normal.iteritems():
    print("normal.{}_{} = {}".format(k, v, v), file=of)

print("\n# set.normal = set", file=of)
for k, v in normal_tumors.iteritems():
    print("set.{} = {} {}".format(k, " ".join(v), k), file=of)
print("\n# set.tumor = set", file=of)
for k, v in tumor_normal.iteritems():
    print("set.{} = {} {}".format(k, " ".join(normal_tumors[v]), v), file=of)
print("\n# set.pair = set", file=of)
for k, v in tumor_normal.iteritems():
    print("set.{}_{} = {} {}".format(k, v, " ".join(normal_tumors[v]), v), file=of)

if args.sample_fastq_file is not None:
    print("\n# sample_fastq_file", file=of)
    sample_fastq = yaml.load(args.sample_fastq_file)
    split_samples = set()
    for k, v in sample_fastq.iteritems():
        for idx, fastq in enumerate(v):
            print("fq.{}.{} = {}".format(k, idx, " ".join(fastq)), file=of)
    for k, v in sample_fastq.iteritems():
        for idx, fastq in enumerate(v):
            print("split.{}.{} = {}".format(k, idx, k), file=of)
            split_samples.add("{}.{}".format(k, idx))
        line = "split.{} = ".format(k)
        for idx, fastq in enumerate(v):
            line += "{}.{} ".format(k, idx)
        print(line.strip(), file=of)
        print("", file=of)
    print("SPLIT_SAMPLES = {}".format(" ".join(split_samples)), file=of)

if args.sample_merge_file is not None:
    print("\n# sample_merge_file", file=of)
    sample_merge = yaml.load(args.sample_merge_file)
    print("MERGE_SAMPLES = {}".format(" ".join(sample_merge.keys())), file=of)
    for k, v in sample_merge.iteritems():
        print("merge.{} = {}".format(k, " ".join(v)), file=of)

print("\n# defaults", file=of)
config = yaml.load(args.config_file)
for k, v in config.iteritems():
    print("{} ?= {}".format(k.upper(), lowerBool(v)), file=of)
