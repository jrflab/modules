#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/../../test_functions.sh
set -e

echo "TEST VARIANT CALLING test_data"
cd ${DIR}/test_data
# create project_config.inc
python modules/configure.py
echo "==somatic_variants: mutect2 + strelka + varscan=="
echo "run"
make -j8 somatic_variants && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/vcf/*mutect*.vcf; do
    compare_vcf $x vcf/$(basename $x)
done
for x in output/vcf/*somatic_indels*.vcf; do
    compare_vcf $x vcf/$(basename $x)
done
echo "==hotspot=="
echo "run"
make -j8 hotspot && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/vcf_ann/*hotspot*.vcf; do
    compare_vcf $x vcf_ann/$(basename $x)
done
echo "==facets copy number=="
echo "run"
make -j8 facets && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/facets/cncf/*.txt; do
    diff_files $x facets/cncf/$(basename $x)
done
for x in output/facets/*.txt; do
    diff_files $x facets/$(basename $x)
done
echo "==somatic annotations=="
echo "run"
make -j8 ann_somatic_vcf && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/vcf_ann/*.vcf; do
    compare_vcf $x vcf_ann/$(basename $x)
done
echo "==mutation summary=="
echo "run"
make -j8 mutation_summary && echo "success" || (echo "failed" && exit 1)
echo "check files"
for x in output/summary/tsv/*.tsv; do
    diff_files $x summary/tsv/$(basename $x)
done
exit 0
