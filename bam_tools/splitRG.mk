include modules/Makefile.inc

LOGDIR = log/splitRG.$(NOW)

split : $(foreach sample,$(SAMPLES),rg/EEC14/$(sample).bam)

define split-rg
rg/EEC14/$1.bam : bam/EEC14.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p rg/EEC14 && \
				      $$(SAMTOOLS) view -b -r $1 bam/EEC14.bam > rg/EEC14/$1.bam && \
				      $$(SAMTOOLS) index rg/EEC14/$1.bam")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call split-rg,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/splitRG.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split
