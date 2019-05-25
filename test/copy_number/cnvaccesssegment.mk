include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_segment.$(NOW)
PHONY += cnvaccess cnvaccess/segmented

cnvaccess_segment : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/segmented/$(sample).RData)

define cnvaccess-segment
cnvaccess/segmented/%.RData : cnvaccess/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p cnvaccess/segmented && \
												 $(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 3 --sample_name $$(*)")
												 
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-segment,$(sample))))
	
.PHONY: $(PHONY)
