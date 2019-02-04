include modules/Makefile.inc

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000'

LOGDIR ?= log/genotype_pdx.$(NOW)
PHONY += sufam summary

MOUSE_SAMPLES ?= $(wildcard $(foreach sample,$(sample_category.mouse),sufam/$(sample).txt))

genotype_pdx : $(foreach sample,$(sample_category.mouse),sufam/$(sample).txt) sufam/PDX.vcf summary/mouse_summary.xlsx summary/tsv/mouse_summary.tsv

sufam/PDX.vcf : summary/tsv/mutation_summary.tsv
	$(call RUN, -c -s 8G -m 16G,"$(RSCRIPT) modules/variant_callers/genotypepdx.R")

define genotype-pdx
sufam/%.txt : bam/%.bam sufam/PDX.vcf
	$$(call RUN,-v $$(SUFAM_ENV) -c -s 2G -m 4G -w 2880,"sufam --sample_name $$(*) $$(SUFAM_OPTS) $$(REF_FASTA) sufam/PDX.vcf bam/$$(*).bam > sufam/$$(*).txt")
	
endef
 $(foreach sample,$(sample_category.mouse),\
		$(eval $(call genotype-pdx,$(sample))))
		
summary/tsv/mouse_summary.tsv summary/mouse_summary.xlsx : $(wildcard $(foreach sample,$(sample_category.mouse),sufam/$(sample).txt))
	$(call RUN,-n 1 -s 4G -m 4G,"$(RSCRIPT) modules/summary/mousesummary.R --in_file '$${sample_category.mouse}' --out_file summary/tsv/mouse_summary.tsv && \
								 python modules/summary/mouse_summary_excel.py")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
