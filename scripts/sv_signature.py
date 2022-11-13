#!/usr/bin/env python

""" extract structural variant signatures
"""

import argparse
import viola
import numpy as np
import pandas as pd
import os
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='sv_signature.py',
                                     description='SV feature extraction')
    parser.add_argument('--bedpe_infile', required=True)
    parser.add_argument('--fragile_bed', required=True)
    parser.add_argument('--timing_bedgraph', required=True)
    parser.add_argument('--sv_definitions', required=True)
    parser.add_argument('--text_outfile', required=True)

    args = parser.parse_args()

    sample_bedpe = viola.viola.read_bedpe(args.bedpe_infile)
    bed_fragile = viola.read_bed(args.fragile_bed)
    bedgraph_timing = viola.read_bed(args.timing_bedgraph)

    sample_bedpe.annotate_bed(bed=bed_fragile, annotation='fragile', how='flag')
    sample_bedpe.annotate_bed(bed=bedgraph_timing, annotation='timing', how='value')
    sample_bedpe.calculate_info('(${timingleft} + ${timingright}) / 2', 'timing')
    
    feature_matrix = sample_bedpe.classify_manual_svtype(definitions=args.sv_definitions, return_data_frame=True)
    feature_matrix.drop('others', axis=1, inplace=True)
    feature_matrix.to_csv(args.text_outfile, index=False, sep='\t')
