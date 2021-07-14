include modules/Makefile.inc

LOGDIR = log/splitRG.$(NOW)

split : $(foreach sample,$(SAMPLES),rg/XXX/$(sample).bam)

define split-rg
rg/XXX/$1.bam : bam/XXX.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p rg/XXX && \
				      $$(SAMTOOLS) view -b -r $1 bam/XXX.bam > rg/XXX/$1.bam && \
				      $$(SAMTOOLS) index rg/XXX/$1.bam")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call split-rg,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/splitRG.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split
