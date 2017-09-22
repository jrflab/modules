#!/usr/bin/env python
""" convert defuse output to unified tsv format
"""

import sys
import argparse
import pandas as pd
import numpy as np

info_fields = '''cluster_id
splitr_sequence
splitr_span_pvalue
splitr_pos_pvalue
splitr_min_pvalue
adjacent
altsplice
break_adj_entropy1
break_adj_entropy2
break_adj_entropy_min
breakpoint_homology
breakseqs_estislands_percident
cdna_breakseqs_percident
deletion
est_breakseqs_percident
eversion
exonboundaries
expression1
expression2
gene_end1
gene_end2
gene_location1
gene_location2
gene_start1
gene_start2
gene_strand1
gene_strand2
genome_breakseqs_percident
genomic_strand1
genomic_strand2
interchromosomal
interrupted_index1
interrupted_index2
inversion
library_name
max_map_count
max_repeat_proportion
mean_map_count
min_map_count
num_multi_map
num_splice_variants
orf
read_through
repeat_proportion1
repeat_proportion2
span_coverage1
span_coverage2
span_coverage_max
span_coverage_min
splice_score
splicing_index1
splicing_index2
upstream_gene
downstream_gene'''.split('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='defuse2usv.py',
                                     description='convert defuse output to unified sv tsv format')
    parser.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'),
                        help='defuse output file')
    args = parser.parse_args()

    df = pd.read_table(args.input)
    df = df.rename(columns={
        'gene_name1': 'GeneSymbol1', 'gene_name2': 'GeneSymbol2',
        'gene1': 'EnsemblGene1', 'gene2': 'EnsemblGene2',
        'gene_chromosome1': 'Chr1', 'gene_chromosome2': 'Chr2',
        'gene_align_strand1': 'Strand1', 'gene_align_strand2': 'Strand2',
        'splitr_count': 'NumSplitReads', 'span_count': 'NumSpanningReads',
        'genomic_break_pos1': 'Breakpoint1', 'genomic_break_pos2': 'Breakpoint2'})

    df['Info'] = ''
    for i, row in df.iterrows():
        fields = []
        for field in info_fields:
            f = row[field]
            if type(f) == str and f == '':
                continue
            if type(f) == float and np.isnan(f):
                continue
            x = "{}=".format(field)
            if ' ' in str(f) or ';' in str(f):
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
