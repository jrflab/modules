#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

OUT_DIR=${DIR}/output
SCRIPT_DIR=${DIR}/../../scripts
TEST_DATA_DIR=${DIR}/test_data
TEST_DATA_OUT_DIR=${DIR}/test_data/output


mkdir -p ${OUT_DIR}
python ${SCRIPT_DIR}/classify_snv_pathogenicity_vcf.py ${TEST_DATA_DIR}/pathogen_test_snps.vcf > ${OUT_DIR}/pathogen_test_out_snps.vcf
compare_files ${OUT_DIR}/pathogen_test_out_snps.vcf ${TEST_DATA_OUT_DIR}/pathogen_test_out_snps.vcf

python ${SCRIPT_DIR}/classify_indel_pathogenicity_vcf.py ${TEST_DATA_DIR}/pathogen_test_indels.vcf > ${OUT_DIR}/pathogen_test_out_indels.vcf
compare_files ${OUT_DIR}/pathogen_test_out_indels.vcf ${TEST_DATA_OUT_DIR}/pathogen_test_out_indels.vcf
