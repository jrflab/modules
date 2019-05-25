include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_plot.$(NOW)
PHONY += cnvkit cnvaccess/plot

cnvaccess_plot : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/plot/$(sample).A.ontarget.pdf
				cnvaccess/log2/$(sample).B.ontarget.pdf
				cnvaccess/log2/$(sample).AB.ontarget.pdf
				cnvaccess/log2/$(sample).offtarget.pdf
				cnvaccess/log2/$(sample).pdf)

define cnvaccess-plot
cnvaccess/log2/%.A.ontarget.pdf cnvaccess/log2/%.B.ontarget.pdf cnvaccess/log2/%.AB.ontarget.pdf cnvaccess/log2/%.offtarget.pdf cnvaccess/log2/%.pdf : cnvaccess/cnr/%.A.cnr cnvaccess/cnr/%.B.cnr cnvaccess/cnr/%.C.cnr
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 1 --sample_name $$(*)")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-plot,$(sample))))
				
.PHONY: $(PHONY)
