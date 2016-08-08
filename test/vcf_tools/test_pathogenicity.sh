#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

OUT_DIR=${DIR}/output
SCRIPT_DIR=${DIR}/../../vcf_tools
TEST_DATA_DIR=${DIR}/test_data
TEST_DATA_OUT_DIR=${DIR}/test_data/output


mkdir -p ${OUT_DIR}
python ${SCRIPT_DIR}/classify_pathogenicity_vcf.py ${TEST_DATA_DIR}/pathogen_test.vcf > ${OUT_DIR}/pathogen_test_out.vcf
compare_files ${OUT_DIR}/pathogen_test_out.vcf ${TEST_DATA_OUT_DIR}/pathogen_test_out.vcf

python ${SCRIPT_DIR}/classify_pathogenicity_vcf.py ${TEST_DATA_DIR}/pathogen_test2.vcf > ${OUT_DIR}/pathogen_test2_out.vcf
compare_files ${OUT_DIR}/pathogen_test2_out.vcf ${TEST_DATA_OUT_DIR}/pathogen_test2_out.vcf
