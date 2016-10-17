#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../../test_functions.sh
set -e

TEST_GENES='ASXL1|ALOX12B|ARID5B|TGFBR1|TGFBR2|TP53|TSC1'

echo "TEST VARIANT CALLING test_data"
cd ${DIR}/test_data
# create project_config.inc
python modules/configure.py
echo "==facets copy number=="
echo "run"
make -j8 facets && echo "success" || (echo "failed" && exit 1)
echo "check files"
# can't check cncf files (non-deterministic)
diff_files <(grep -P "$TEST_GENES" output/facets/geneCN.txt) <(grep -P "$TEST_GENES" facets/geneCN.txt)
