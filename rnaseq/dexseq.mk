include modules/Makefile.inc

LOGDIR ?= log/exon_counts.$(NOW)
PHONY += dex_seq

dex_seq : $(foreach sample,$(TUMOR_SAMPLES),dex_seq/$(sample).taskcomplete)

define exon-count
dex_seq/%.txt : star/bam/%.star.sorted.filtered.bam
	$$(call RUN,-c -s 8G -m 12G,"echo $$<")
dex_seq/%.taskcomplete : dex_seq/%.txt
	$$(call RUN,-c -s 1G -m 1G,"echo $$<")
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call exon-count,$sample)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
