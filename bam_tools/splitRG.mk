include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

split : $(foreach sample,$(SAMPLES),bam/EEC3/$(sample).bam) \
	$(foreach sample,$(SAMPLES),bam/EEC3/$(sample).bam.bai)

define split-rg
bam/EEC3/$1.bam : etc/bam/EEC3-2.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p bam/EEC3 && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@)")

bam/EEC3/$1.bam.bai : bam/EEC3/$1.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call split-rg,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/split_rg.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split
