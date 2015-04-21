# annotate titan module
include modules/Makefile.inc
include modules/copy_number/titan.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/ann_titan.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : ann_all ann_varscan ann_mutect ann_strelka

ann_all : ann_varscan ann_mutect ann_strelka

ann_varscan : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach type,varscan_indels varscan_snps,\
	$(foreach suff,$(call VCF_SUFFIXES,$(type)),vcf/$(pair).$(suff).titan.vcf)))

ann_strelka : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach type,strelka_indels strelka_snps,\
	$(foreach suff,$(call VCF_SUFFIXES,$(type)),vcf/$(pair).$(suff).titan.vcf)))

ann_mutect : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call VCF_SUFFIXES,mutect),vcf/$(pair).$(suff).titan.vcf))

vcf/%.titan.vcf : vcf/%.vcf titan/optclust_results_w$(TITAN_WINDOW_SIZE)_p$(DEFAULT_PLOIDY_PRIOR)/%.titan_seg.txt
	$(call LSCRIPT_MEM,4G,6G,"$(ANNOTATE_TITAN_LOH_VCF) --titanSeg $(<<) --outFile $@ $<")
