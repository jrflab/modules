include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_plot_test.$(NOW)
PHONY += cnvkit cnvaccess/plot

cnvkit : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/plot/$(sample).A.ontarget.pdf cnvaccess/plot/$(sample).B.ontarget.pdf cnvaccess/plot/$(sample).AB.ontarget.pdf cnvaccess/plot/$(sample).offtarget.pdf)

define cnvkit-plot
cnvaccess/plot/%.A.ontarget.pdf cnvaccess/plot/%.B.ontarget.pdf cnvaccess/plot/%.AB.ontarget.pdf cnvaccess/plot/%.offtarget.pdf : cnvaccess/cnr/%.timestamp
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/cnvkitplot.R --tumor $$(*) --normals '$$(NORMAL_SAMPLES)'")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-plot,$(sample))))
				
.PHONY: $(PHONY)
