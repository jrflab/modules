include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_plot.$(NOW)
PHONY += cnvkit cnvaccess/plot cnvaccess/plot/log2

cnvaccess_plot : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/plot/log2/$(sample).pdf)

define cnvaccess-plot
cnvaccess/plot/log2/%.pdf : cnvaccess/cnr/%.A.cnr cnvaccess/cnr/%.B.cnr cnvaccess/cnr/%.C.cnr
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/cnvaccess.R --type 1 --sample_name $$(*)")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-plot,$(sample))))
				
.PHONY: $(PHONY)
