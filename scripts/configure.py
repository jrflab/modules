#!/usr/bin/env python

from __future__ import print_function
import yaml
import argparse
import collections

def lowerBool(x):
    if isinstance(x, bool):
        return str(x).lower()
    else:
        return x

def sample_yaml2mk(samples_file, out):
    samples = yaml.full_load(open(args.samples_file, 'r'))

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

    print("", file=out)
    print("SAMPLES = {} {}".format(" ".join(tumors), " ".join(normals)), file=out)
    print("TUMOR_SAMPLES = {}".format(" ".join(tumors)), file=out)
    print("NORMAL_SAMPLES = {}".format(" ".join(normals)), file=out)
    print("SAMPLE_PAIRS = ", end="", file=out)
    for k, v in tumor_normal.items():
        print("{}_{}".format(k, v), end=" ", file=out)
    print("", file=out)
    print("SAMPLE_SETS = ", end="", file=out)
    for normal, tumors in normal_tumors.items():
        print("{tumors}_{normal}".format(tumors="_".join(tumors), normal=normal), end=" ", file=out)
    print("", file=out)

    print("\n# tumor.normal = tumors", file=out)
    for normal, tumors in normal_tumors.items():
        print("tumor.{normal} = {tumors}".format(normal=normal, tumors=" ".join(tumors)), file=out)
    print("\n# normal.tumor = normal", file=out)
    for tumor, normal in tumor_normal.items():
        print("normal.{tumor} = {normal}".format(tumor=tumor, normal=normal), file=out)
    print("\n# tumor.pair = tumor", file=out)
    for tumor, normal in tumor_normal.items():
        print("tumor.{tumor}_{normal} = {tumor}".format(tumor=tumor, normal=normal), file=out)
    print("\n# tumors.set = tumors", file=out)
    for normal, tumors in normal_tumors.items():
        print("tumors.{tumors}_{normal} = {tumors_space}".format(tumors="_".join(tumors), normal=normal,
                                                                 tumors_space=" ".join(tumors)), file=out)
    print("\n# normal.pair = normal", file=out)
    for tumor, normal in tumor_normal.items():
        print("normal.{tumor}_{normal} = {normal}".format(tumor=tumor, normal=normal), file=out)
    print("\n# normal.set = normal", file=out)
    for normal, tumors in normal_tumors.items():
        print("normal.{tumors}_{normal} = {normal}".format(tumors="_".join(tumors), normal=normal), file=out)

    print("\n# set.normal = set", file=out)
    for normal, tumors in normal_tumors.items():
        print("set.{normal} = {tumors} {normal}".format(normal=normal, tumors=" ".join(tumors)), file=out)
    print("\n# set.tumor = set", file=out)
    for tumor, normal in tumor_normal.items():
        print("set.{tumor} = {tumors} {normal}".format(tumor=tumor, tumors=" ".join(normal_tumors[normal]),
                                                       normal=normal), file=out)
    print("\n# set.pair = set", file=out)
    for tumor, normal in tumor_normal.items():
        print("set.{tumor}_{normal} = {set}".format(tumor=tumor, normal=normal,
                                                    set=" ".join(normal_tumors[normal])), file=out)

    print("\n# set.set = set", file=out)
    for normal, tumors in normal_tumors.items():
        print("set.{tumors}_{normal} = {set}".format(tumors="_".join(tumors), normal=normal,
                                                       set=" ".join(tumors + [normal])),
              file=out)

    print("\n# pairs.set = pairs", file=out)
    for normal, tumors in normal_tumors.items():
        print("pairs.{tumors}_{normal} = {pairs}".format(tumors="_".join(tumors), normal=normal,
                                                         pairs=" ".join(
                                                             ["_".join([tumor, normal]) for tumor in tumors])),
              file=out)


def sample_attr_yaml2mk(sample_attr_file, out):
    print("\n# sample_attr_file", file=out)
    sample_attr = yaml.full_load(open(sample_attr_file, 'r'))
    for attr, m in sample_attr.items():
        for k, v in m.items():
            print("{}.{} = {}".format(attr, k, v), file=out)


def sample_fastq_yaml2mk(sample_fastq_file, out):
    print("\n# sample_fastq_file", file=out)
    sample_fastq = yaml.full_load(open(sample_fastq_file, 'r'))
    split_samples = set()
    for k, v in sample_fastq.items():
        for idx, fastq in enumerate(v):
            print("fq.{}.{} = {}".format(k, idx, " ".join(fastq)), file=out)
    for k, v in sample_fastq.items():
        for idx, fastq in enumerate(v):
            print("split.{}.{} = {}".format(k, idx, k), file=out)
            split_samples.add("{}.{}".format(k, idx))
        line = "split.{} = ".format(k)
        for idx, fastq in enumerate(v):
            line += "{}.{} ".format(k, idx)
        print(line.strip(), file=out)
        print("", file=out)
    print("SPLIT_SAMPLES = {}".format(" ".join(split_samples)), file=out)


def sample_merge_yaml2mk(sample_merge_file, out):
    print("\n# sample_merge_file", file=out)
    sample_merge = yaml.full_load(open(args.sample_merge_file, 'r'))
    print("MERGE_SAMPLES = {}".format(" ".join(list(sample_merge.keys()))), file=out)
    for k, v in sample_merge.items():
        print("merge.{} = {}".format(k, " ".join(v)), file=out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='configure', description='Convert project YAML file to Make')
    parser.add_argument('--project_config_file', help='project yaml config file', default='project_config.yaml')
    parser.add_argument('--samples_file', help='yaml samples file', default='samples.yaml')
    parser.add_argument('--sample_attr_file', help='yaml sample attr file', default='sample_attr.yaml')
    parser.add_argument('--sample_fastq_file', help='yaml sample fastq file mappings', default='sample.fastq.yaml')
    parser.add_argument('--sample_merge_file', help='yaml sample merge mappings')
    parser.add_argument('--out_file', help='project make include file', nargs='?', default='project_config.inc')
    args = parser.parse_args()

    of = open(args.out_file, 'w')

    config = yaml.full_load(open(args.project_config_file, 'r'))
    for k, v in config.items():
        print("{} = {}".format(k.upper(), lowerBool(v)), file=of)

    try:
        sample_yaml2mk(args.samples_file, of)
    except IOError:
        print("Error loading {}, skipping".format(args.samples_file))

    try:
        sample_attr_yaml2mk(args.sample_attr_file, of)
    except:
        print("Error loading {}, skipping".format(args.sample_attr_file))

    if args.sample_fastq_file is not None:
        try:
            sample_fastq_yaml2mk(args.sample_fastq_file, of)
        except IOError:
            print("Error loading {}, skipping.".format(args.sample_fastq_file))

    if args.sample_merge_file is not None:
        try:
            sample_merge_yaml2mk(args.sample_merge_file, of)
        except IOError:
            print("Error loading {}, skipping.".format(args.sample_merge_file))

    of.close()
