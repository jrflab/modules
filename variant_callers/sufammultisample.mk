include modules/Makefile.inc

LOGDIR ?= log/sufam_multisample.$(NOW)
PHONY += sufam summary

SUFAM_SUMMARY ?= $(wildcard $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv))

sufam_multisample : $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv) summary/sufam_summary.xlsx

ifeq ($(GATK_SPLIT_CHR),true

define combine-samples
sufam/%.txt : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/combinesamples.R --sample_set $$*")

sufam/%.tsv : sufam/%.txt
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/updatesamples.R --sample_set $$*")
	
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples,$(set))))

summary/sufam_summary.xlsx : $(wildcard $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv))
	$(call RUN,-c -s 12G -m 16G,"export R_LIBS='/lila/data/reis-filho/usr/anaconda-envs/jrflab-modules-0.1.4/lib/R/library:/lila/data/reis-filho/usr/lib64/R/library' &&\
								 $(RSCRIPT) modules/summary/sufamsummary.R --in_file '$(SUFAM_SUMMARY)'")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)