include modules/Makefile.inc

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000'

LOGDIR ?= log/genotype_pdx.$(NOW)
PHONY += sufam

summary_to_vcf : sufam/PDX.vcf
genotype_pdx : $(foreach sample,$(sample_category.mouse),sufam/$(sample).txt)

sufam/PDX.vcf : summary/tsv/mutation_summary.tsv
	$(call RUN, -c -s 8G -m 16G,"$(RSCRIPT) modules/variant_callers/genotypepdx.R")

define genotype-pdx
sufam/%.txt : bam/%.bam sufam/PDX.vcf
	$$(call RUN,-v $$(SUFAM_ENV) -c -s 2G -m 4G -w 2880,"sufam --sample_name $$(*) $$(SUFAM_OPTS) $$(REF_FASTA) sufam/PDX.vcf bam/$$(*).bam > sufam/$$(*).txt")
	
endef
 $(foreach sample,$(sample_category.mouse),\
		$(eval $(call genotype-pdx,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

