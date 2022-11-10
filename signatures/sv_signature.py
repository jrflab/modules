#!/usr/bin/env python

""" compute structural variant signatures
"""

import argparse
import viola
import numpy as np
import pandas as pd
import os
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='sv_signature.py', description='Compute structural variant signatures')
    parser.add_argument('--pcawg_bedpe', nargs='?', default='/data/reis-filho/lib/resource_files/viola/pcawg/', help='Path containing PCAWG bedpe files')
    parser.add_argument('--fragile_site', nargs='?', default='/data/reis-filho/lib/resource_files/viola/annotation/fragile_site.hg19.bed', help='Fragile sites')
    parser.add_argument('--replication_timing', nargs='?', default='/data/reis-filho/lib/resource_files/viola/annotation/replication_timing.bedgraph', help='Replication timing')
    parser.add_argument('--sv_definition', nargs='?', default='/data/reis-filho/lib/resource_files/viola/definitions/sv_class_definition.txt', help='SV definition')
    parser.add_argument('--feature_matrix', nargs='?', default='feature_maxtrix.txt', help='Feature matrix path')
    parser.add_argument('--exposure_matrix', nargs='?', default='exposure_maxtrix.txt', help='Exposure matrix path')
    parser.add_argument('--signature_matrix', nargs='?', default='signature_maxtrix.txt', help='Signature matrix path')
    parser.add_argument('--name', nargs='?', default='viola-sv', help='Run name')

    args = parser.parse_args()

    pcawg_bedpe = viola.read_bedpe_multi(args.pcawg_bedpe)
    bed_fragile = viola.read_bed(args.fragile_site)
    bedgraph_timing = viola.read_bed(args.replication_timing)
    pcawg_bedpe.annotate_bed(bed=bed_fragile, annotation='fragile', how='flag')
    pcawg_bedpe.annotate_bed(bed=bedgraph_timing, annotation='timing', how='value')
    pcawg_bedpe.calculate_info('(${timingleft} + ${timingright}) / 2', 'timing')
    feature_matrix = pcawg_bedpe.classify_manual_svtype(definitions=args.pcawg_bedpe, return_data_frame=True)
    feature_matrix.drop('others', axis=1, inplace=True)
    result_silhouette, result_metrics, exposure_matrix, signature_matrix = viola.SV_signature_extractor(feature_matrix, n_iter=10, name=args.name, n_components=12, init='nndsvda', solver='mu', beta_loss='kullback-leibler', max_iter=10000, random_state=1)

    feature_matrix.to_csv(args.feature_matrix, index=False, sep='\t')
    
    fh = open(args.exposure_matrix, 'r+')
    np.savetxt(fh, exposure_matrix)
    fh.close()

    fh = open(args.signature_matrix, 'r+')
    np.savetxt(fh, signature_matrix)
    fh.close()