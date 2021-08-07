include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

split : $(foreach sample,$(SAMPLES),bam/HEC6-ISHI/$(sample).bam) \
	$(foreach sample,$(SAMPLES),bam/HEC6-ISHI/$(sample).bam.bai)

define split-rg
bam/HEC6-ISHI/$1.bam : etc/bam/HEC6-ISHI-2.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p bam/HEC6-ISHI && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@)")

bam/HEC6-ISHI/$1.bam.bai : bam/HEC6-ISHI/$1.bam
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
