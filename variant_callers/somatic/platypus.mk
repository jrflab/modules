include modules/Makefile.inc
include modules/variant_callers/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)

PLATYPUS_ENV = $(HOME)/share/usr/anaconda-envs/platypus-0.8.1

PHONY += platypus_indels
platypus_indels: $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).platypus_indels.vcf)

define platypus-tumor-normal
vcf/$1_$2.platypus_indels.vcf : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_ENV_PARALLEL_MEM,$$(PLATYPUS_ENV),"platypus callVariants --bamFiles=$$(<)$$(,)$$(<<) --nCPU 4 --refFile=$$(REF_FASTA) --output=$$@")
endef

include modules/vcf_tools/vcftools.mk
