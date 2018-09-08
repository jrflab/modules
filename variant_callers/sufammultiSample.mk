include modules/Makefile.inc

LOGDIR ?= log/sufam_multisample.$(NOW)
PHONY += sufam

sufal_multisample : $(foreach sample,$(NORMAL_SAMPLES),sufam/$(sample).tsv)

define combine-samples
sufam/%.tsv : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/combineSamples.R --patient $$*")
	
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call combine-samples,$(sample))))
