include modules/Makefile.inc

LOGDIR = log/splitRG.$(NOW)

split : $(foreach sample,$(SAMPLES),rg/EEC128/$(sample).bam)

define split-rg
rg/EEC128/$1.bam : bam/EEC128.bam
	$$(call RUN,-n 1 -s 4G -m 6G,"set -o pipefail && \
				      mkdir -p rg/EEC128 && \
				      $$(SAMTOOLS) view -b -r $1 bam/EEC128.bam > rg/EEC128/$1.bam")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-fastq,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/splitRG.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split
