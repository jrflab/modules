# annotate titan module
include modules/Makefile.inc
include modules/copy_number/titan.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/ann_titan.$(NOW)

#VCF_FIELDS += titanCN titanMinorCN titanMajorCN titanCall titanMedianRatio titanMedianLogR

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : ann_all ann_varscan ann_mutect ann_strelka

ann_all : ann_varscan ann_mutect ann_strelka

ann_varscan : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach type,varscan_indels varscan_snps,\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$(type)),vcf/$(pair).$(suff).titan.vcf))) \
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach type,varscan_indels varscan_snps,\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$(type)),\
	$(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),\
	tables/$(pair).$(suff).titan.$(ext).txt)))) \
$(foreach type,varscan_indels varscan_snps,\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$(type)),\
	$(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),alltables/allTN.$(suff).titan.$(ext).txt)))

ann_strelka : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach type,strelka_indels strelka_snps,\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$(type)),vcf/$(pair).$(suff).titan.vcf))) \
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach type,strelka_indels strelka_snps,\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$(type)),\
	$(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),\
	tables/$(pair).$(suff).titan.$(ext).txt)))) \
$(foreach type,strelka_indels strelka_snps,\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$(type)),\
	$(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),alltables/allTN.$(suff).titan.$(ext).txt)))

ann_mutect : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,mutect),vcf/$(pair).$(suff).titan.vcf)) \
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,mutect),\
	$(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),\
	tables/$(pair).$(suff).titan.$(ext).txt))) \
$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,mutect),\
	$(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),alltables/allTN.$(suff).titan.$(ext).txt))

$(info vcf/$1.%.titan.vcf : vcf/$1.%.vcf titan/optclust_results_w$(TITAN_WINDOW_SIZE)_p$(DEFAULT_PLOIDY_PRIOR)/$1.titan_seg.txt)
define titan-pair
vcf/$1.%.titan.vcf : vcf/$1.%.vcf titan/optclust_results_w$(TITAN_WINDOW_SIZE)_p$(DEFAULT_PLOIDY_PRIOR)/$1.titan_seg.txt
	$$(call LSCRIPT_MEM,4G,6G,"$$(ANNOTATE_TITAN_LOH_VCF) --titanSeg $$(<<) --outFile $$@ $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call titan-pair,$(pair))))

include modules/vcf_tools/vcftools.mk
