include modules/Makefile.inc

LOGDIR ?= log/sufam_multisample.$(NOW)
PHONY += sufam summary

SUFAM_SUMMARY ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),sufam/$(sample).tsv))

sufam_multisample : $(foreach sample,$(NORMAL_SAMPLES),sufam/$(sample).tsv) summary/sufam_summary.xlsx

define combine-samples
sufam/%.txt : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/combineSamples.R --patient $$*")

sufam/%.tsv : sufam/%.txt
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/updateSamples.R --patient $$*")
	
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call combine-samples,$(sample))))

summary/sufam_summary.xlsx : $(wildcard sufam/$(NORMAL_SAMPLES).tsv)
	$(call RUN,-c -s 12G -m 16G,"export R_LIBS='/lila/data/reis-filho/usr/anaconda-envs/jrflab-modules-0.1.4/lib/R/library:/lila/data/reis-filho/usr/lib64/R/library' &&\
								 $(RSCRIPT) modules/summary/sufamsummary.R --in_file '$(SUFAM_SUMMARY)'")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)