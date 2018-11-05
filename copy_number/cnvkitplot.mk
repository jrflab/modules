include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_plot.$(NOW)
PHONY += cnvkit cnvkit/plot

cnvkit : $(foreach sample,$(SAMPLES),cnvkit/plot/$(sample).ontarget.pdf cnvkit/plot/$(sample).offtarget.pdf)

define cnvkit-plot
cnvkit/plot/%.ontarget.pdf cnvkit/plot/%.offtarget.pdf : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 4G -m 6G,"$(RSCRIPT) modules/copy_number/cnvkitplot.R --in_file $$(<)")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvkit-plot,$(sample))))
				
.PHONY: $(PHONY)
