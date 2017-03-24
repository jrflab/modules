# Run multiple indel callers and then merge them
#
include modules/Makefile.inc
LOGDIR = log/somatic_variants.$(NOW)

SNV_TYPE ?= mutect
INDEL_TYPE ?= strelka_varscan_indels
VARIANT_TYPES ?= $(SNV_TYPE) $(INDEL_TYPE)

PHONY += all_somatic
all_somatic: somatic_vcfs somatic_tables facets

CONCAT_VCF = python modules/vcf_tools/concat_vcf.py

vcf/%.somatic_indels.vcf : vcf/%.mutect_indels.vcf vcf/%.strelka_varscan_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"$(CONCAT_VCF) $^ | $(VCF_SORT) $(REF_DICT) - > $@")

vcf/%.strelka_varscan_indels.vcf : vcf/%.varscan_indels.vcf vcf/%.strelka_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"(grep -P '^#' $<; $(BEDTOOLS) intersect -a $< -b <($(PASS_FILTER_VCF) $(<<))) | uniq > $@")

include modules/variant_callers/somatic/mutect.mk
include modules/variant_callers/somatic/mutect2.mk
include modules/variant_callers/somatic/strelka.mk
include modules/variant_callers/somatic/varscanTN.mk
include modules/copy_number/facets.mk
include modules/vcf_tools/annotateSomaticVcf.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY) 

