#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

OUT_DIR=${DIR}/output
SCRIPT_DIR=${DIR}/../../scripts
TEST_DATA_DIR=${DIR}/test_data
TEST_DATA_OUT_DIR=${DIR}/test_data/output

mkdir -p ${OUT_DIR}
python ${SCRIPT_DIR}/create_sample_yaml.py --sample_file ${OUT_DIR}/samples.yaml --sample_fastq_file ${OUT_DIR}/sample.fastq.yaml ${TEST_DATA_DIR}/rawdata
compare_files ${OUT_DIR}/samples.yaml ${TEST_DATA_OUT_DIR}/samples.yaml
#compare_files ${OUT_DIR}/sample.fastq.yaml ${TEST_DATA_OUT_DIR}/sample.fastq.yaml will fail on travis CI bc of abs paths
