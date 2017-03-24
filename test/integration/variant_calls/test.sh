#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../../test_functions.sh
set -e

echo "TEST VARIANT CALLING test_data"
cd ${DIR}/test_data
# create project_config.inc
python modules/configure.py
echo "==somatic_variants: mutect + strelka + varscan + facets + annotations=="
echo "run"
make -j8 somatic_variants && echo "success" || (echo "failed" && exit 1)
echo "==hotspot=="
echo "run"
make -j8 hotspot && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/vcf/*.vcf; do
    compare_vcf $x vcf/$(basename $x)
done
for x in output/vcf_ann/*.vcf; do
    compare_vcf $x vcf_ann/$(basename $x)
done
echo "==mutation summary=="
echo "run"
make -j8 mutation_summary && echo "success" || (echo "failed" && exit 1)
#echo "check files"
#for x in output/summary/tsv/*.tsv; do
    #diff_files $x summary/tsv/$(basename $x)
#done
exit 0
