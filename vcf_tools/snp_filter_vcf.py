#!/usr/bin/env python
""" output snps only
"""

import vcf
import sys

if __name__ == "__main__":
    vcf_reader = vcf.Reader(sys.stdin)
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    for record in vcf_reader:
        if record.is_snp:
            vcf_writer.write_record(record)

    vcf_writer.close()
