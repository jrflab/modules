include modules/Makefile.inc
include modules/clonality/setuppyclone.mk

LOGDIR ?= log/plot_pyclone.$(NOW)
PHONY += pyclone

plot_pyclone : $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/report/pyclone.pdf)

define plot-pyclone
pyclone/%/report/pyclone.pdf : pyclone/%/report/pyclone.tsv
	$$(call RUN,-s 4G -m 6G -w 7200,"$(RSCRIPT) modules/clonality/plotpyclone.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) --min_depth $(MIN_DEPTH)")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call plot-pyclone,$(set))))
