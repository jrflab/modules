include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_segment.$(NOW)
PHONY += cnvaccess cnvaccess/segmented cnvaccess/plot/segmented

cnvaccess_segment : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/segmented/$(sample).RData cnvaccess/plot/segmented/$(sample).pdf)


define cnvaccess-segment
cnvaccess/segmented/%.RData : cnvaccess/cnr/%.pool-A.cnr cnvaccess/cnr/%.pool-B.cnr cnvaccess/cnr/%.no-pool.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p cnvaccess/segmented && \
												 $(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 3 --sample_name $$(*)")
												 
cnvaccess/plot/segmented/%.pdf : cnvaccess/segmented/%.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p cnvaccess/plot/segmented && \
												 $(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 4 --sample_name $$(*)")
												 
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-segment,$(sample))))
	
.PHONY: $(PHONY)
