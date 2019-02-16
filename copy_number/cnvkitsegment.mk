include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_segment.$(NOW)
PHONY += cnvkit cnvkit/segmented cnvkit/totalcopy cnvkit/called

cnvkit_segment : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/totalcopy/$(sample).RData) $(foreach sample,$(TUMOR_SAMPLES),cnvkit/segmented/$(sample).pdf) $(foreach sample,$(TUMOR_SAMPLES),cnvkit/called/$(sample).RData)

define cnvkit-totalcopy
cnvkit/segmented/%.pdf cnvkit/totalcopy/%.RData : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p cnvkit/segmented && \
												 mkdir -p cnvkit/totalcopy && \
												 $(RSCRIPT) modules/copy_number/cnvkit.R --type total-copy --sample_name $$(*)")
												 
cnvkit/called/%.RData : cnvkit/totalcopy/%.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p cnvkit/called && \
												 $(RSCRIPT) modules/copy_number/cnvkit.R --type call-cna --sample_name $$(*)")

endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-totalcopy,$(sample))))
	
.PHONY: $(PHONY)
