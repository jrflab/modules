#!/usr/bin/env python
from __future__ import print_function
import yaml
import argparse
import collections

""" convert yaml files to make include files
"""


def lowerBool(x):
    if isinstance(x, bool):
        return str(x).lower()
    else:
        return x


def sample_yaml2mk(samples_file, out_file):
    samples = yaml.load(open(args.samples_file, 'r'))

    tumors = set()
    normals = set()
    normal_tumors = collections.defaultdict(list)
    tumor_normal = {}
    for s in samples:
        if 'normal' in s:
            normals.add(s['normal'])
        if 'tumor' in s:
            for t in s['tumor']:
                tumors.add(t)
        if 'normal' in s and 'tumor' in s:
            normal_tumors[s['normal']].extend(s['tumor'])
            for t in s['tumor']:
                tumor_normal[t] = s['normal']

    print("", file=out_file)
    print("SAMPLES = {} {}".format(" ".join(tumors), " ".join(normals)), file=out_file)
    print("TUMOR_SAMPLES = {}".format(" ".join(tumors)), file=out_file)
    print("NORMAL_SAMPLES = {}".format(" ".join(normals)), file=out_file)
    print("SAMPLE_PAIRS = ", end="", file=out_file)
    for k, v in tumor_normal.iteritems():
        print("{}_{}".format(k, v), end=" ", file=out_file)
    print("", file=out_file)
    print("SAMPLE_SETS = ", end="", file=out_file)
    for k, v in normal_tumors.iteritems():
        print("{}_{}".format("_".join(v), k), end=" ", file=out_file)
    print("", file=out_file)

    print("\n# tumor.normal = tumors", file=out_file)
    for k, v in normal_tumors.iteritems():
        print("tumor.{} = {}".format(k, " ".join(v)), file=out_file)
    print("\n# normal.tumor = normal", file=out_file)
    for k, v in tumor_normal.iteritems():
        print("normal.{} = {}".format(k, v), file=out_file)
    print("\n# tumor.pair = tumor", file=out_file)
    for k, v in tumor_normal.iteritems():
        print("tumor.{}_{} = {}".format(k, v, k), file=out_file)
    print("\n# normal.pair = normal", file=out_file)
    for k, v in tumor_normal.iteritems():
        print("normal.{}_{} = {}".format(k, v, v), file=out_file)

    print("\n# set.normal = set", file=out_file)
    for k, v in normal_tumors.iteritems():
        print("set.{} = {} {}".format(k, " ".join(v), k), file=out_file)
    print("\n# set.tumor = set", file=out_file)
    for k, v in tumor_normal.iteritems():
        print("set.{} = {} {}".format(k, " ".join(normal_tumors[v]), v), file=out_file)
    print("\n# set.pair = set", file=out_file)
    for k, v in tumor_normal.iteritems():
        print("set.{}_{} = {} {}".format(k, v, " ".join(normal_tumors[v]), v), file=out_file)


def sample_fastq_yaml2mk(sample_fastq_file, out_file):
    print("\n# sample_fastq_file", file=out_file)
    sample_fastq = yaml.load(open(sample_fastq_file, 'r'))
    split_samples = set()
    for k, v in sample_fastq.iteritems():
        for idx, fastq in enumerate(v):
            print("fq.{}.{} = {}".format(k, idx, " ".join(fastq)), file=out_file)
    for k, v in sample_fastq.iteritems():
        for idx, fastq in enumerate(v):
            print("split.{}.{} = {}".format(k, idx, k), file=out_file)
            split_samples.add("{}.{}".format(k, idx))
        line = "split.{} = ".format(k)
        for idx, fastq in enumerate(v):
            line += "{}.{} ".format(k, idx)
        print(line.strip(), file=out_file)
        print("", file=out_file)
    print("SPLIT_SAMPLES = {}".format(" ".join(split_samples)), file=out_file)


def sample_merge_yaml2mk(sample_merge_file, out_file):
    print("\n# sample_merge_file", file=out_file)
    sample_merge = yaml.load(open(args.sample_merge_file, 'r'))
    print("MERGE_SAMPLES = {}".format(" ".join(sample_merge.keys())), file=out_file)
    for k, v in sample_merge.iteritems():
        print("merge.{} = {}".format(k, " ".join(v)), file=out_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='configure',
                                     description='Convert project YAML file to Make')
    parser.add_argument('--project_config_file', help='project yaml config file',
                        default='project_config.yaml')
    parser.add_argument('--samples_file', help='yaml samples file', default='samples.yaml')
    parser.add_argument('--sample_fastq_file', help='yaml sample fastq file mappings', default='sample.fastq.yaml')
    parser.add_argument('--sample_merge_file', help='yaml sample merge mappings')
    parser.add_argument('--out_file', help='project make include file', type=argparse.FileType('w', 0), nargs='?',
                        default='project_config.inc')
    args = parser.parse_args()

    of = args.out_file

    config = yaml.load(open(args.project_config_file, 'r'))
    for k, v in config.iteritems():
        print("{} = {}".format(k.upper(), lowerBool(v)), file=of)

    try:
        sample_yaml2mk(args.samples_file, args.out_file)
    except IOError:
        print("Error loading {}, skipping".format(args.samples_file))

    if args.sample_fastq_file is not None:
        try:
            sample_fastq_yaml2mk(args.sample_fastq_file, args.out_file)
        except IOError:
            print("Error loading {}, skipping.".format(args.sample_fastq_file))

    if args.sample_merge_file is not None:
        try:
            sample_merge_yaml2mk(args.sample_merge_file, args.out_file)
        except IOError:
            print("Error loading {}, skipping.".format(args.sample_merge_file))
