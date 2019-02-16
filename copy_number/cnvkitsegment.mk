include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_segment.$(NOW)
PHONY += cnvkit cnvkit/segmented

cnvkit_segment : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/segmented/$(sample).RData)

define cnvkit-segment
cnvkit/segmented/%.RData : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"$(RSCRIPT) modules/copy_number/cnvkit.R --type segment --sample_name $$(*)")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-segment,$(sample))))
				
.PHONY: $(PHONY)

