# merge split fastqs for workflows like defuse

include modules/Makefile.inc

LOGDIR ?= log/merge_split_fastq.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : fastq

fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz fastq/$(sample).2.fastq.gz)

define merged-fastq
fastq/$1.%.fastq.gz : $$(foreach split,$2,fastq/$$(split).%.fastq.gz)
	$$(call LSCRIPT,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-fastq,$(sample),$(split_lookup.$(sample)))))

