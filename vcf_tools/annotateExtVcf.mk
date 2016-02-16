include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/vcf.$(NOW)

EXT_NAME ?= ext

SOMATIC_FILTERS := 

ext_ann : ext_vcfs ext_tables

ext_vcfs : $(call SOMATIC_VCFS,$(EXT_NAME))
ext_tables : $(call SOMATIC_TABLES,$(EXT_NAME))

include modules/vcf_tools/vcftools.mk
