#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../../test_functions.sh
set -e

echo "TEST VARIANT CALLING test_data"
cd ${DIR}/test_data
# create project_config.inc
python modules/configure.py
echo "==mutect=="
echo "run"
make mutect && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/alltables/*mutect*; do
    compare_files $x alltables/$(basename $x)
done
echo "==strelka=="
echo "run"
make strelka && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/alltables/*.strelka_{snps,indels}.*; do
    compare_files $x alltables/$(basename $x)
done
echo "==varscanTN=="
echo "run"
make varscanTN && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/alltables/*.varscan_{snps,indels}.*; do
    compare_files $x alltables/$(basename $x)
done
echo "==merge_strelka_varscan=="
echo "run"
make merge_strelka_varscan && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/alltables/*strelka_varscan*; do
    compare_files $x alltables/$(basename $x)
done
echo "==hotspotTN=="
echo "run"
make hotspotTN && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/alltables/*hotspot*; do
    compare_files $x alltables/$(basename $x)
done
make mutation_summary && echo "success" || (echo "failed" && exit 1)
echo "check files"
compare_files summary/mutation_summary.xlsx output/summary/mutation_summary.xlsx
exit 0
