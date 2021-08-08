include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

split : $(foreach sample,$(SAMPLES),bam/EEC91/$(sample).bam) \
	$(foreach sample,$(SAMPLES),bam/EEC91/$(sample).bam.bai)

define split-rg
bam/EEC91/$1.bam : etc/bam/EEC91-1.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p bam/EEC91 && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@)")

bam/EEC91/$1.bam.bai : bam/EEC91/$1.bam
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
