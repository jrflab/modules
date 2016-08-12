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
for x in output/vcf/*mutect*.vcf; do
    compare_vcf $x vcf/$(basename $x)
done
echo "==strelka_varscan_indels=="
echo "run"
make strelka_varscan_indels && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/vcf/*.strelka_varscan_indels.vcf; do
    compare_vcf $x vcf/$(basename $x)
done
echo "==hotspot=="
echo "run"
make hotspot && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/vcf_ann/*hotspot*.vcf; do
    compare_vcf $x vcf_ann/$(basename $x)
done
make ann_somatic_vcf && echo "success" || (echo "failed" && exit 1)
for x in output/vcf_ann/*.vcf; do
    compare_vcf $x vcf_ann/$(basename $x)
done
make mutation_summary && echo "success" || (echo "failed" && exit 1)
echo "check files"
compare_files output/summary/mutation_summary.xlsx summary/mutation_summary.xlsx
exit 0
