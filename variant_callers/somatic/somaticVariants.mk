include modules/Makefile.inc
LOGDIR = log/somatic_variants.$(NOW)

SNV_TYPE ?= mutect
INDEL_TYPE ?= somatic_indels
VARIANT_TYPES ?= $(SNV_TYPE) $(INDEL_TYPE)

PHONY += all_somatic
all_somatic: somatic_vcfs somatic_tables facets

CONCAT_VCF = python modules/vcf_tools/concat_vcf.py

vcf/%.somatic_variants.vcf : vcf/%.$(SNV_TYPE).vcf vcf/%.$(INDEL_TYPE).vcf
	$(call RUN,-s 9G -m 12G,"$(CONCAT_VCF) $^ | $(VCF_SORT) $(REF_DICT) - > $@")

include modules/variant_callers/somatic/mutect.mk
include modules/variant_callers/somatic/somaticIndels.mk
include modules/copy_number/facets.mk
include modules/vcf_tools/annotateSomaticVcf.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY) 

