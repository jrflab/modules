include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

split : $(foreach sample,$(SAMPLES),bam/EEC98/$(sample).bam) \
	$(foreach sample,$(SAMPLES),bam/EEC98/$(sample).bam.bai)

define split-rg
bam/EEC98/$1.bam : etc/bam/EEC98.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p bam/EEC98 && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@)")

bam/EEC98/$1.bam.bai : bam/EEC98/$1.bam
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
