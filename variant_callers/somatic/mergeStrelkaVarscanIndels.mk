# Merge strelka and varscan indel results
LOGDIR ?= log/merge_strelka_varscan_indels.$(NOW)

include modules/Makefile.inc
include modules/config.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : strelka_varscan_merge strelka_varscan_merge_vcfs 

strelka_varscan_merge : strelka_varscan_merge_vcfs 
strelka_varscan_merge_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).strelka_varscan_indels.vcf)

vcf_ann/%.strelka_varscan_indels.vcf : vcf_ann/%.varscan_indels.vcf vcf_ann/%.strelka_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"grep -P '^#' $< > $@ && $(BEDTOOLS) intersect -a $< -b $(<<) >> $@")

include modules/vcf_tools/vcftools.mk
