# Run multiple indel callers and then merge them
#
include modules/Makefile.inc
LOGDIR = log/strelkaVarscan.$(NOW)


.PHONY : merge_indel_vcfs
strelka_varscan_merge_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).somatic_indels.vcf)
#strelka_varscan_merge_mafs : $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).strelka_varscan_indels.vcf)

CONCAT_VCF = python modules/vcf_tools/concat_vcf.py

vcf/%.somatic_indels.vcf : vcf/%.mutect_indels.vcf vcf/%.strelka_varscan_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"$(CONCAT_VCF) $^ | $(VCF_SORT) $(REF_DICT) - > $@")

vcf/%.strelka_varscan_indels.vcf : vcf/%.varscan_indels.vcf vcf/%.strelka_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"grep -P '^#' $< > $@ && $(BEDTOOLS) intersect -a $< -b $(<<) >> $@")


include modules/variant_callers/somatic/strelka.mk
include modules/variant_callers/somatic/varscanTN.mk
