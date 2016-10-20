#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../../test_functions.sh
set -e

TEST_GENES='\tSYCP2\t|'

echo "TEST VARIANT CALLING test_data"
cd ${DIR}/test_data
# create project_config.inc
python modules/configure.py
echo "==facets copy number=="
echo "run"
make -j8 facets && echo "success" || (echo "failed" && exit 1)
echo "check files"
# can't check cncf files (non-deterministic)
for x in output/facets/cncf/*.cncf.txt; do
    compare_cncf $x facets/cncf/$(basename $x)
done
