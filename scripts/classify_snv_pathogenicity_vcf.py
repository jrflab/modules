#!/usr/bin/env python

""" classify pathogenicity of vcf records, querying provean as necessary
"""

import argparse
import vcf
import sys
import classify_pathogenicity_vcf as cp


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='classify_pathogenicity_vcf.py',
                                     description='Add pathogenicity to vcf file')
    parser.add_argument('vcf_infile')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    assert "ANN" in vcf_reader.infos
    assert "HOTSPOT" in vcf_reader.infos or "hotspot" in vcf_reader.infos
    assert "fathmm_pred" in vcf_reader.infos
    assert "facetsLOH" in vcf_reader.infos
    assert "MutationTaster_pred" in vcf_reader.infos

    # add necessary info headers
    vcf_reader.infos['pathogenicity'] = vcf.parser._Info(id='pathogenicity', num=-1, type='String',
                                                         desc="Classification of pathogenicity",
                                                         source=None, version=None)
    records = [x for x in vcf_reader]
    for record in records:
        cp.classify_pathogenicity(record)

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    for record in records:
        vcf_writer.write_record(record)
    vcf_writer.close()
