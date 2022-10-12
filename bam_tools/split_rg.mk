include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

split_rg : $(foreach sample,$(SAMPLES),split_rg/$(sample).bam) \
	   $(foreach sample,$(SAMPLES),split_rg/$(sample).bam.bai)

define split-rg
split_rg/$1.bam : bam/
	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@)")

split_rg/$1.bam.bai : split_rg/$1.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<) && \
				      cp $$(@) split_rg/$$(*).bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call split-rg,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/split_rg.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: split_rg
