include modules/Makefile.inc
include modules/variant_callers/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)

PLATYPUS_ENV = $(HOME)/share/usr/anaconda-envs/platypus-0.8.1

PHONY += platypus_indels
platypus_indels: $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).platypus_indels.vcf)

define platypus-tumor-normal-chr
platypus/chr_vcf/$1_$2.$3.platypus.vcf : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_ENV_PARALLEL_MEM,$$(PLATYPUS_ENV),4,2G,3G,"platypus callVariants --regions=$3 \
		--bamFiles=$$(<)$$(,)$$(<<) --nCPU 4 --refFile=$$(REF_FASTA) --output=$$@")
endef
$(foreach chr,$(CHROMOSOMES),\
	$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call platypus-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

INDEL_FILTER_VCF = python modules/vcf_tools/indel_filter_vcf.py
SNP_FILTER_VCF = python modules/vcf_tools/snp_filter_vcf.py

vcf/%.platypus_indels.vcf : $(foreach chr,$(CHROMOSOMES),platypus/chr_vcf/%.$(chr).platypus.vcf)
	$(call LSCRIPT_MEM,4G,8G,"(grep '^#' $<; cat $^ | grep -v '^#' | \
		$(VCF_SORT) $(REF_DICT) - ) | $(INDEL_FILTER_VCF) > $@")

vcf/%.platypus_snps.vcf : $(foreach chr,$(CHROMOSOMES),platypus/chr_vcf/%.$(chr).platypus.vcf)
	$(call LSCRIPT_MEM,4G,8G,"(grep '^#' $<; cat $^ | grep -v '^#' | \
		$(VCF_SORT) $(REF_DICT) - ) | $(SNP_FILTER_VCF) > $@")

include modules/vcf_tools/vcftools.mk
