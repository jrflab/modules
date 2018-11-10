include modules/Makefile.inc

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000'

LOGDIR ?= log/genotype_hotspots.$(NOW)
PHONY += hotspot

genotype_hotspots : $(foreach sample,$(SAMPLES),hotspot/$(sample).txt)

define genotype-hotspots
hotspot/%.txt : bam/%.bam
	$$(call RUN,-v $$(SUFAM_ENV) -c -s 2G -m 4G -w 2880,"sufam --sample_name $$(*) $$(SUFAM_OPTS) $$(REF_FASTA) modules/reference/hotspots/hotspot-dedup.vcf bam/$$(*).bam > hotspot/$$(*).txt")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-hotspots,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

