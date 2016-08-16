include modules/Makefile.inc

LOGDIR ?= log/vcf.$(NOW)

EXT_NAME ?= ext

SOMATIC_FILTERS := 

ext_ann : ext_vcfs ext_tables

ext_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).ext.vcf)
ext_tables : altables/allTN.ext.tab.txt

include modules/vcf_tools/vcftools.mk
