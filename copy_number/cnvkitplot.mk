include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_plot.$(NOW)
PHONY += cnvkit

cnvkit : $(foreach sample,$(SAMPLES),cnvkit/$(sample).ontarget.pdf cnvkit/$(sample).offtarget.pdf)

define cnvkit-plot
cnvkit/%.ontarget.pdf cnvkit/%.offtarget.pdf : cnvkit/%.cnr
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/copy_number/cnvkitplot.R --in_file $$(<)")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvkit-plot,$(sample))))
				
.PHONY: $(PHONY)
