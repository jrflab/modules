include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_plot_test.$(NOW)
PHONY += cnvkit cnvkit/plot

cnvkit : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/plot/$(sample).A.ontarget.pdf cnvkit/plot/$(sample).B.ontarget.pdf cnvkit/plot/$(sample).AB.ontarget.pdf cnvkit/plot/$(sample).offtarget.pdf)

define cnvkit-plot
cnvkit/plot/%.A.ontarget.pdf cnvkit/plot/%.B.ontarget.pdf cnvkit/plot/%.AB.ontarget.pdf cnvkit/plot/%.offtarget.pdf : cnvkit/cnr/%.timestamp
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/cnvkitplot.R --tumor $$(*) --normals '$$(NORMAL_SAMPLES)'")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-plot,$(sample))))
				
.PHONY: $(PHONY)
