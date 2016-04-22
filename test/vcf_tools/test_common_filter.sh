#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

mkdir -p ${DIR}/output
python ${DIR/test\//}/common_filter_vcf.py ${DIR}/test_data/test.vcf ${DIR}/output/testout.vcf
compare_files ${DIR}/output/testout.vcf ${DIR}/test_data/output/testout.vcf
