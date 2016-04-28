# Run VarScan and strelka on tumour-normal matched pairs for indels
#
include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc
LOGDIR = log/strelkaVarscan.$(NOW)


.PHONY : strelka_varscan_merge_vcfs strelka_varscan_merge
strelka_varscan_merge : strelka_varscan_merge_vcfs 
strelka_varscan_merge_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).strelka_varscan_indels.vcf)
#strelka_varscan_merge_mafs : $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).strelka_varscan_indels.vcf)

vcf_ann/%.strelka_varscan_indels.vcf : vcf_ann/%.varscan_indels.vcf vcf_ann/%.strelka_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"grep -P '^#' $< > $@ && $(BEDTOOLS) intersect -a $< -b $(<<) >> $@")


include modules/variant_callers/somatic/strelka.mk
include modules/variant_callers/somatic/varscanTN.mk
