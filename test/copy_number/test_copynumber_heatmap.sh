#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../test_functions.sh
set -e

mkdir -p ${DIR}/output
Rscript ${DIR/test\//}/facetsGeneCNPlot.R ${DIR}/test_data/geneCN.txt ${DIR}/output/copynumber_heatmap.pdf
compare_files ${DIR}/output/copynumber_heatmap.pdf ${DIR}/test_data/output/copynumber_heatmap.pdf
