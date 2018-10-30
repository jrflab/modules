include modules/Makefile.inc

DEFAULT_ENV = /home/$(USER)/share/usr/anaconda-envs/jrflab-modules-0.1.5

LOGDIR ?= log/sufam_multisample.$(NOW)
PHONY += sufam

sufam_multisample : $(foreach sample,$(NORMAL_SAMPLES),sufam/$(sample).tsv)

define combine-samples
sufam/%.txt : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/combineSamples.R --patient $$*")

sufam/%.tsv : sufam/%.txt
	$$(call RUN,-c -s 4G -m 6G --env $$(DEFAULT_ENV),"$(RSCRIPT) modules/variant_callers/updateSamples.R --patient $$*")
	
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call combine-samples,$(sample))))
