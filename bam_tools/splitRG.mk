include modules/Makefile.inc

LOGDIR = log/splitRG.$(NOW)

split : $(foreach sample,$(SAMPLES),rg/EEC87/$(sample).bam)

define split-rg
rg/EEC87/$1.bam : bam/EEC87.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p rg/EEC87 && \
				      $$(SAMTOOLS) view -b -r $1 bam/EEC87.bam > rg/EEC87/$1.bam && \
				      $$(SAMTOOLS) index rg/EEC87/$1.bam")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call split-rg,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/splitRG.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split
