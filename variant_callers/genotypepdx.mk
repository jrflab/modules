include modules/Makefile.inc

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000'

LOGDIR ?= log/genotype_pdx.$(NOW)
PHONY += sufam

genotype_pdx : $(foreach sample,$(sample_category.mouse),sufam/$(sample).txt)

define genotype-pdx
sufam/%.txt : bam/%.bam
	$$(call RUN,-v $$(SUFAM_ENV) -c -s 2G -m 4G -w 2880,"sufam --sample_name $$(*) $$(SUFAM_OPTS) $$(REF_FASTA) modules/reference/hotspots/hotspot-dedup.vcf bam/$$(*).bam > sufam/$$(*).txt")
	
endef
 $(foreach sample,$(sample_category.mouse),\
		$(eval $(call genotype-pdx,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

