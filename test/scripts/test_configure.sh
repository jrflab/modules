#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

OUT_DIR=${DIR}/output
SCRIPT_DIR=${DIR}/../../
TEST_DATA_DIR=${DIR}/test_data
TEST_DATA_OUT_DIR=${DIR}/test_data/output

mkdir -p ${DIR}/output
python ${SCRIPT_DIR}/configure.py --project_config_file ${DIR}/../../default_project_config.yaml --samples_file ${TEST_DATA_DIR}/samples.yaml --sample_fastq_file ${TEST_DATA_DIR}/sample.fastq.yaml --out_file ${OUT_DIR}/project_config.inc
compare_files ${OUT_DIR}/project_config.inc ${TEST_DATA_OUT_DIR}/project_config.inc

