include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

split : $(foreach sample,$(SAMPLES),bam/EEC131/$(sample).bam) \
	$(foreach sample,$(SAMPLES),bam/EEC131/$(sample).bam.bai)

define split-rg
bam/EEC131/$1.bam : etc/bam/EEC131.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p bam/EEC131 && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@)")

bam/EEC131/$1.bam.bai : bam/EEC131/$1.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call split-rg,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/split_rg.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split
