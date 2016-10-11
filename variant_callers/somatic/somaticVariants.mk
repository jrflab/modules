# Run multiple indel callers and then merge them
#
include modules/Makefile.inc
LOGDIR = log/somatic_variants.$(NOW)

SNV_TYPE ?= mutect_snps
INDEL_TYPE ?= somatic_indels

.PHONY : somatic_variants somatic_indels somatic_snvs
somatic_variants: somatic_snvs somatic_indels
somatic_snvs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(SNV_TYPE).vcf)
somatic_indels : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(INDEL_TYPE).vcf)
#strelka_varscan_merge_mafs : $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).strelka_varscan_indels.vcf)

CONCAT_VCF = python modules/vcf_tools/concat_vcf.py

vcf/%.somatic_indels.vcf : vcf/%.mutect_indels.vcf vcf/%.strelka_varscan_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"$(CONCAT_VCF) $^ | $(VCF_SORT) $(REF_DICT) - > $@")

vcf/%.strelka_varscan_indels.vcf : vcf/%.varscan_indels.vcf vcf/%.strelka_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"(grep -P '^#' $<; $(BEDTOOLS) intersect -a $< -b $(<<)) | uniq > $@")

include modules/variant_callers/somatic/mutect.mk
include modules/variant_callers/somatic/mutect2.mk
include modules/variant_callers/somatic/strelka.mk
include modules/variant_callers/somatic/varscanTN.mk
