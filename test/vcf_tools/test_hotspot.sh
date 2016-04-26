#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

mkdir -p ${DIR}/output
python ${DIR/test\//}/hotspot_vcf.py ${DIR}/test_data/hotspot.vcf \
    ${DIR}/test_data/cancer_hotspots.txt \
    > ${DIR}/output/hotspot.vcf
compare_files ${DIR}/output/hotspot.vcf ${DIR}/test_data/output/hotspot.vcf
