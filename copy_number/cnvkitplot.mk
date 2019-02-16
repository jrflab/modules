include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_plot.$(NOW)
PHONY += cnvkit cnvkit/log2

cnvkit_plot : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/log2/$(sample).ontarget.pdf cnvkit/log2/$(sample).offtarget.pdf)

define cnvkit-plot
cnvkit/log2/%.ontarget.pdf cnvkit/log2/%.offtarget.pdf : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 4G -m 6G,"$(RSCRIPT) modules/copy_number/cnvkitplot.R --in_file $$(<)")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-plot,$(sample))))
				
.PHONY: $(PHONY)
