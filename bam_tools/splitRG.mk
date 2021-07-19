include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

split : $(foreach sample,$(SAMPLES),bam/EEC25/$(sample).bam)

define split-rg
bam/EEC25/$1.bam : etc/bam/EEC25-1.bam
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      mkdir -p bam/EEC25 && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@) && \
				      $$(SAMTOOLS) index $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call split-rg,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/split_rg.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split
