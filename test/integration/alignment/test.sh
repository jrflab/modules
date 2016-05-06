#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

echo "TEST ALIGNMENT TARGETED SEQUENCING test_data"
cd ${DIR}/test_data
# create project_config.inc
python modules/configure.py
echo "==bwamem=="
echo "run"
make bwamem && echo "success" || (echo "failed" && exit 1)
echo "check files"
compare_bams bam/test1T.bam output/bam/test1T.bam

exit 0
