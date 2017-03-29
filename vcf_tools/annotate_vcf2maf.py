#!/usr/bin/env python
# run vcf2maf on a vcf file and then merge in the annotations

import argparse
import vcf
import sys
import pandas as pd
import tempfile
import subprocess
import os
import re


def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='annotate_vcf2maf.py',
                                     description='annotate a vcf file with vcf2maf maf annotations')
    parser.add_argument('--vcf2maf', default='/opt/common/CentOS_6-dev/perl/perl-5.22.0/bin/perl '
                        '/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/vcf2maf.pl',
                        help='vcf2maf perl command')
    parser.add_argument('--vcf2maf_opts', default='--vep-data /opt/common/CentOS_6-dev/vep/v86/ --ncbi-build GRCh37'
                        ' --maf-center mskcc.org --tmp-dir /scratch/tmpbvHWSR'
                        ' --vep-path /opt/common/CentOS_6-dev/vep/v86/'
                        ' --custom-enst /opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_at_mskcc'
                        ' --species homo_sapiens', help='vcf2maf options')
    parser.add_argument('--ref_fasta', default='/ifs/depot/assemblies/H.sapiens/b37_dmp/b37.fasta',
                        help='reference fasta')
    parser.add_argument('--vep_forks', default=1, type=int, help='number of VEP forks')
    parser.add_argument('--filter_vcf',
                        default='/opt/common/CentOS_6-dev/vep/v86/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz',
                        help='ExAC non-TCGA vcf file')
    parser.add_argument('--hotspot_list',
                        default=None,
                        help='hotspot list file')
    parser.add_argument('vcf_infile')
    args = parser.parse_args()

    # tidy hotspot list
    hotspots = None
    if args.hotspot_list is not None:
        hdf = pd.read_table(args.hotspot_list, sep='\t')
        hdf = tidy_split(hdf, 'Variants')
        hdf['Variant'] = hdf['Variants'].str.split(':').apply(lambda x: x[0])
        hgvsp = []
        for i, row in hdf.iterrows():
            if re.match(r'^[A-Z]', row['Residue']):
                hgvsp.append('p.' + row['Residue'] + row['Variant'])
            else:
                hgvsp.append('p.' + row['Variant'])
        hdf['HGVSp'] = hgvsp
        hotspots = set(hdf['Gene'] + ":" + hdf['HGVSp'])

    maf_file = '{}/{}'.format(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))

    # run vcf2maf
    FNULL = open(os.devnull, 'w')
    cmd = '{vcf2maf} {vcf2maf_opts}' \
        ' --ref-fasta {ref_fasta}' \
        ' --input-vcf {vcf}' \
        ' --filter-vcf {filter_vcf}' \
        ' --vep-forks {vep_forks}' \
        ' --output-maf {maf}'.format(vcf2maf=args.vcf2maf, vcf2maf_opts=args.vcf2maf_opts,
                                     ref_fasta=args.ref_fasta, vcf=args.vcf_infile,
                                     vep_forks=args.vep_forks, filter_vcf=args.filter_vcf,
                                     maf=maf_file)
    retcode = subprocess.call(cmd, shell=True, stdout=FNULL, stderr=FNULL)
    if retcode != 0:
        sys.stderr.write(cmd + '\n')
        sys.stderr.write('vcf2maf failed\n')
        sys.exit(1)

    maf = pd.read_table(maf_file, dtype={'Chromosome': str}, comment='#')
    os.remove(maf_file)

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))
    vcf_reader.infos['HGVSp_Short'] = vcf.parser._Info(id='HGVSp_Short', num='.',
                                                       type='String',
                                                       desc='HGVSp short from vcf2maf',
                                                       source='vcf2maf',
                                                       version=None)
    vcf_reader.infos['SYMBOL'] = vcf.parser._Info(id='SYMBOL', num='.',
                                                  type='String',
                                                  desc='gene symbol from vcf2maf',
                                                  source='vcf2maf',
                                                  version=None)
    vcf_reader.infos['Variant_Classification'] = vcf.parser._Info(id='Variant_Classification', num='.',
                                                                  type='String',
                                                                  desc='variant classification from vcf2maf',
                                                                  source='vcf2maf',
                                                                  version=None)

    if hotspots is not None:
        vcf_reader.infos['cmo_hotspot'] = vcf.parser._Info(id='cmo_hotspot', num='0',
                                                           type='Flag',
                                                           desc='cmo hotspot',
                                                           source='vcf2maf',
                                                           version=None)

    records = [x for x in vcf_reader]
    for i, record in enumerate(records):
        assert record.CHROM == maf.ix[i, 'Chromosome']
        record.INFO['HGVSp_Short'] = maf.ix[i, 'HGVSp_Short']
        record.INFO['SYMBOL'] = maf.ix[i, 'SYMBOL']
        record.INFO['Variant_Classification'] = maf.ix[i, 'Variant_Classification']
        gene_hgvsp = str(record.INFO['SYMBOL']) + ':' + str(record.INFO['HGVSp_Short'])
        if hotspots is not None and gene_hgvsp in hotspots:
            record.INFO['cmo_hotspot'] = True
    assert i == len(maf) - 1

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    for record in records:
        vcf_writer.write_record(record)
    vcf_writer.close()
