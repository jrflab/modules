include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_plot.$(NOW)
PHONY += cnvkit cnvaccess/plot cnvaccess/plot/log2 cnvaccess/plot/segmented cnvaccess/plot/bychr

cnvaccess_plot : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/plot/log2/$(sample).pdf cnvaccess/plot/bychr/$(sample)/timestamp)

define cnvaccess-plot
cnvaccess/plot/log2/%.pdf : cnvaccess/cnr/%.pool-A.cnr cnvaccess/cnr/%.pool-B.cnr cnvaccess/cnr/%.no-pool.cnr
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 1 --sample_name $$(*)")
	
cnvaccess/plot/bychr/%/timestamp : cnvaccess/cnr/%.pool-A.cnr cnvaccess/cnr/%.pool-B.cnr cnvaccess/cnr/%.no-pool.cnr
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"mkdir -p cnvaccess/plot/ && \
																	 mkdir -p cnvaccess/plot/bychr/ && \
																	 mkdir -p cnvaccess/plot/bychr/$$(*)/ && \
																	 $(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 2 --sample_name $$(*)")

endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-plot,$(sample))))
				
.PHONY: $(PHONY)
