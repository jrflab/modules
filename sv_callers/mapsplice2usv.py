#!/usr/bin/env python
""" convert mapsplice output to unified tsv format
"""

import sys
import argparse
import pandas as pd
import numpy as np
import gffutils
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='mapsplice2usv.py',
                                     description='convert mapsplice output to unified sv tsv format')
    parser.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'),
                        help='mapsplice output file')
    parser.add_argument('--genes_db', required=True, help='genes db file')
    args = parser.parse_args()

    db = gffutils.FeatureDB(args.genes_db, keep_order=True)

    df = pd.read_table(args.input, header=None)
    df = df.drop(62, axis=1)
    cols = ["Chr12", "Breakpoint1", "Breakpoint2", "id", "coverage", "strand", "rgb", "block_count", "block_size",
            "block_distance", "entropy", "flank_case", "flank_string", "min_mismatch", "max_mismatch", "ave_mismatch",
            "max_min_suffix", "max_min_prefix", "min_anchor_difference", "unique_read_count", "multi_read_count",
            "paired_read_count", "left_paired_read_count", "right_paired_read_count", "multiple_paired_read_count",
            "unique_paired_read_count", "single_read_count", "encompassing_read_pair_count", "donor_start",
            "acceptor_end", "donor_isoforms", "acceptor_isoforms", "obsolete1", "obsolete2", "obsolete3",
            "obsolete4", "minimal_donor_isoform_length", "maximal_donor_isoform_length",
            "minimal_acceptor_isoform_length", "maximal_acceptor_isoform_length", "paired_reads_entropy",
            "mismatch_per_bp", "anchor_score", "max_donor_fragment", "max_acceptor_fragment",
            "max_cur_fragment", "min_cur_fragment", "ave_cur_fragment", "donor_encompass_unique",
            "donor_encompass_multiple", "acceptor_encompass_unique", "acceptor_encompass_multiple",
            "donor_match_to_normal", "acceptor_match_to_normal", "donor_seq", "acceptor_seq",
            "match_gene_strand", "annotated_type", "fusion_type", "gene_strand", "annotated_gene_donor",
            "annotated_gene_acceptor"]
    df.columns = cols
    df['Chr1'], df['Chr2'] = df['Chr12'].str.split('~', 1).str
    df['Strand1'] = df['strand'].str[0]
    df['Strand2'] = df['strand'].str[1]
    df['GeneSymbol1'] = None
    df['GeneSymbol2'] = None
    df['EnsemblGene1'] = None
    df['EnsemblGene2'] = None
    df = df.rename(columns={'encompassing_read_pair_count': 'NumSpanningReads', 'multi_read_count': 'NumSplitReads'})

    for i, row in df.iterrows():
        gene1_ids = set()
        gene1_symbols = set()
        gene2_ids = set()
        gene2_symbols = set()
        features = list(db.region(seqid=row['Chr1'], start=row['Breakpoint1'],
                                  end=row['Breakpoint1'], strand=row['Strand1'],
                                  completely_within=False, featuretype = 'gene'))
        for feature in features:
            gene1_ids = gene1_ids.union(feature.attributes['gene_id'])
            gene1_symbols = gene1_symbols.union(feature.attributes['gene_name'])
        features = list(db.region(seqid=row['Chr2'], start=row['Breakpoint2'],
                                  end=row['Breakpoint2'], strand=row['Strand2'],
                                  completely_within=False, featuretype = 'gene'))
        for feature in features:
            gene2_ids = gene2_ids.union(feature.attributes['gene_id'])
            gene2_symbols = gene2_symbols.union(feature.attributes['gene_name'])
        df.ix[i, 'GeneSymbol1'] = '|'.join(gene1_symbols)
        df.ix[i, 'GeneSymbol2'] = '|'.join(gene2_symbols)
        df.ix[i, 'EnsemblGene1'] = '|'.join(gene1_ids)
        df.ix[i, 'EnsemblGene2'] = '|'.join(gene2_ids)

    info_fields = ["id", "coverage", "rgb", "block_count", "block_size",
                   "block_distance", "entropy", "flank_case", "flank_string", "min_mismatch", "max_mismatch",
                   "ave_mismatch", "max_min_suffix", "max_min_prefix", "min_anchor_difference",
                   "unique_read_count",
                   "paired_read_count", "left_paired_read_count", "right_paired_read_count", "multiple_paired_read_count",
                   "unique_paired_read_count", "single_read_count", "donor_start",
                   "acceptor_end", "donor_isoforms", "acceptor_isoforms",
                   "minimal_donor_isoform_length", "maximal_donor_isoform_length",
                   "minimal_acceptor_isoform_length", "maximal_acceptor_isoform_length", "paired_reads_entropy",
                   "mismatch_per_bp", "anchor_score", "max_donor_fragment", "max_acceptor_fragment",
                   "max_cur_fragment", "min_cur_fragment", "ave_cur_fragment", "donor_encompass_unique",
                   "donor_encompass_multiple", "acceptor_encompass_unique", "acceptor_encompass_multiple",
                   "donor_match_to_normal", "acceptor_match_to_normal", "donor_seq", "acceptor_seq"]
    df['Info'] = ''
    for i, row in df.iterrows():
        fields = []
        for field in info_fields:
            if type(row[field]) != str or row[field] != '':
                x = "{}=".format(field)
                f = re.sub(r'[\|,]$', r'', str(row[field]))
                if ' ' in str(f):
                    x += '"{}"'.format(f)
                else:
                    x += '{}'.format(f)
                fields.append(x)
        df.ix[i, 'Info'] = ";".join(fields)
    usv_cols = ['Chr1', 'Breakpoint1', 'Strand1', 'Chr2', 'Breakpoint2', 'Strand2',
                'NumSplitReads', 'NumSpanningReads',
                'GeneSymbol1', 'GeneSymbol2',
                'EnsemblGene1', 'EnsemblGene2', 'Info']
    usv_df = df[usv_cols].copy()
    usv_df.to_csv(sys.stdout, index=False, sep='\t')
