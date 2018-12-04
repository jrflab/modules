include modules/Makefile.inc

LOGDIR ?= log/mask_regions.$(NOW)
PHONY += masked

mask : $(foreach sample,$(SAMPLES),masked/$(sample).bam)

define bedtools-mask
masked/%.bam : bam/%.bam
	$$(call RUN,-c -s 6G -m 8G,"bedtools intersect -abam $$(<) -b $$(ONTARGET_FILE) -v > mask/$$(*).bam && \
								samtools index mask/$$(*).bam && \
								cp mask/$$(*).bai mask/$$(*).bam.bai")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call bedtools-mask,$(sample))))
		
.PHONY: $(PHONY)
