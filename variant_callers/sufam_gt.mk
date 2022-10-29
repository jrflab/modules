include modules/Makefile.inc

LOGDIR ?= log/sufam_gt.$(NOW)

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000'

sufam_gt : $(foreach set,$(SAMPLE_SETS),sufam/$(set).vcf)

define sufam-genotype
sufam/$1.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					 --option 1 \
					 --sample_set $1 \
					 --normal_sample '$(NORMAL_SAMPLES)' \
					 --input_file $$(<) \
					 --output_file $$(@)")



#hotspot/%.txt : bam/%.bam
#	$$(call RUN,-v $$(SUFAM_ENV) -c -s 2G -m 4G -w 2880,"sufam --sample_name $$(*) $$(SUFAM_OPTS) $$(REF_FASTA) modules/reference/hotspots/hotspot-dedup.vcf bam/$$(*).bam > hotspot/$$(*).txt")
	
endef
 $(foreach set,$(SAMPLE_SETS),\
		$(eval $(call sufam-genotype,$(set))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

