#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

OUT_DIR=${DIR}/output
SCRIPT_DIR=${DIR}/../../vcf_tools
TEST_DATA_DIR=${DIR}/test_data
TEST_DATA_OUT_DIR=${DIR}/test_data/output

mkdir -p ${OUT_DIR}
python ${SCRIPT_DIR}/merge_vcf.py --out_file ${OUT_DIR}/META31T_META31N.varscan_snps.ft.vcf ${TEST_DATA_DIR}/META31T_META31N.varscan_snps.som_ad_ft.vcf ${TEST_DATA_DIR}/META31T_META31N.varscan_snps.target_ft.vcf
diff_files ${OUT_DIR}/META31T_META31N.varscan_snps.ft.vcf ${TEST_DATA_OUT_DIR}/META31T_META31N.varscan_snps.ft.vcf
python ${SCRIPT_DIR}/merge_vcf.py --out_file ${OUT_DIR}/META31T_META31N.varscan_indels.ann2.vcf ${TEST_DATA_DIR}/META31T_META31N.varscan_indels.ft2.mut_taste.vcf
diff_files ${OUT_DIR}/META31T_META31N.varscan_indels.ann2.vcf ${TEST_DATA_OUT_DIR}/META31T_META31N.varscan_indels.ann2.vcf
