#!/usr/bin/env python
"""
Add VCF INFO field of type FLAG that indicates whether the mutation is a known
hotspot mutation. Requires annotated VCF file, using ANN column. Output new vcf
to stdout.
"""
import argparse
import vcf
import pandas as pd
import sys


def one_to_three_amino_acid_code(x):
    # 3 > 1 AA code mapping
    d = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
     'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}
    # inverse to get 1 > 3 AA code mapping
    inv_d = {v: k for k, v in list(d.items())}
    if x == "*":
        return "*"
    else:
        return inv_d[x]


def read_hotspot(tsv):
    """Read hotspot tsv, group by Hugo Symbol and make aa_change_list attribute
    with three letter code amino acid list."""
    def canonicalize_variant_amino_acids(variants):
        try:
            # string looks like A:5|S:17|A:21, where numbers are counts
            # some _splice annotations have SS:19, so only use len(AA) == 1
            rv = [one_to_three_amino_acid_code(aa_count.split(":")[0]) for
                  aa_count in variants.split("|")
                  if len(aa_count.split(":")[0]) == 1]
        except AttributeError:
            # variant column is NA for _splice values
            rv = []
        return rv

    df = pd.read_csv(tsv, sep="\t")
    df["aa_change_list"] = df[["Codon", "Variant Amino Acid"]]\
        .apply(lambda x:
               ["p.{}{}{}".format(one_to_three_amino_acid_code(x.Codon[0]),
                                  x.Codon[1:], aa) for aa in
                canonicalize_variant_amino_acids(x["Variant Amino Acid"])],
               axis=1)
    # flatten aa changes to single list http://stackoverflow.com/questions/952914
    aa_changes_per_gene = df.groupby("Hugo Symbol").aa_change_list.apply(lambda l: [item for sublist in l for item in sublist])
    return aa_changes_per_gene


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_infile')
    parser.add_argument('hotspot_tsv', help="hotspot file from cancerhotspots.org")

    args = parser.parse_args()

    # parse hotspot file
    aa_changes_per_gene = read_hotspot(args.hotspot_tsv)

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    # Get ANN.Gene and ANN.HGVS.p column indices
    ann_columns = "".join(vcf_reader.infos["ANN"].desc.split(":")[1:]).split("|")
    gene_index, hgvsp_index = [i for i, x in enumerate(ann_columns) if x.strip() == "Gene_Name" or x.strip() == "HGVS.p"]
    assert(ann_columns[gene_index].strip() == "Gene_Name")
    assert(ann_columns[hgvsp_index].strip() == "HGVS.p")

    # add hotspot info fields
    # boolean flag
    vcf_reader.infos['hotspot'] = vcf.parser._Info(id='hotspot', num=0,
                                                   type='Flag',
                                                   desc='Hotspot in cancerhotspots.org',
                                                   source='http://cancerhotspots.org/',
                                                   version=None)
    # per annotation
    vcf_reader.infos['hotspot_per_annotation'] = vcf.parser._Info(id='hotspot', num='.',
                                                   type='String',
                                                   desc='Hotspot in cancerhotspots.org',
                                                   source='http://cancerhotspots.org/',
                                                   version=None)

    # add hotspot info fields to outputted records
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    for record in vcf_reader:
        is_hotspot = False
        hotspots = []

        for ann in record.INFO['ANN']:
            gene = ann.split('|')[gene_index]
            hgvsp = ann.split('|')[hgvsp_index]
            if gene in aa_changes_per_gene and \
               hgvsp in aa_changes_per_gene[gene]:
                hotspots += ['1']
                is_hotspot = True
            else:
                hotspots += ['0']

        record.INFO['hotspot_per_annotation'] = hotspots
        record.INFO['hotspot'] = is_hotspot

        vcf_writer.write_record(record)
    vcf_writer.close()
