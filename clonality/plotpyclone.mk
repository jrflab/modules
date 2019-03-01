include modules/Makefile.inc

LOGDIR ?= log/plot_pyclone.$(NOW)
PHONY += pyclone

plot_pyclone : $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/plots/2_by_2_scatter_plots.pdf)

define plot-pyclone
pyclone/%/plots/2_by_2_scatter_plots.pdf : pyclone/%/pyclone.tsv
	$$(call RUN,-s 4G -m 6G -w 7200,"$(RSCRIPT) modules/clonality/plotpyclone.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) --burnin 5000")

endef
$(foreach set,$(SAMPEL_SETS),\
		$(eval $(call plot-pyclone,$(set))))
