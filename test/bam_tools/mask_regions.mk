include modules/Makefile.inc

LOGDIR ?= log/mask_regions.$(NOW)
PHONY += masked

mask : $(foreach sample,$(SAMPLES),masked/$(sample).bam)

define bedtools-mask
masked/%.bam : bam/%.bam
	$$(call RUN,-c -s 6G -m 8G -w 7200,"bedtools intersect -abam $$(<) -b $$(ONTARGET_FILE) -v > masked/$$(*).bam && \
								samtools index masked/$$(*).bam && \
								cp masked/$$(*).bam.bai masked/$$(*).bai")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call bedtools-mask,$(sample))))
		
.PHONY: $(PHONY)
