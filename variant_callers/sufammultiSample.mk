include modules/Makefile.inc

LOGDIR ?= log/sufal_multisample.$(NOW)
PHONY += sufam_multisample

sufal_multisample : $(foreach sample,$(NORMAL_SAMPLES),sufam_multisample/$(sample).tsv)

define combine-samples
sufam_multisample/%.tsv : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/combineSamples.R --patient $$<")
	
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call combine-samples,$(sample))))
