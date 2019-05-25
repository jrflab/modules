include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_segment.$(NOW)
PHONY += cnvaccess cnvaccess/segmented

cnvaccess_segment : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/segmented/$(sample).RData)

define cnvkit-segment
cnvkit/segmented/%.RData : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p cnvkit/segmented && \
												 $(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 3 --sample_name $$(*)")
												 
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-segment,$(sample))))
	
.PHONY: $(PHONY)
