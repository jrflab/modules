include modules/Makefile.inc

LOGDIR ?= log/sufam_multisample.$(NOW)
PHONY += sufam summary

sufam_multisample : $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv) summary/sufam_summary.xlsx

ifeq ($(PDX),true)

define combine-samples-pdx
sufam/%.txt : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/combinesamplesf.R --sample_set $$*")

sufam/%.tsv : sufam/%.txt
	$$(call RUN,-c -s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/updatesamples.R --sample_set $$*")
	
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples-pdx,$(set))))
		
else 

define combine-samples
sufam/%.txt : summary/tsv/mutation_summary.tsv
	$$(call RUN,-s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/combinesamples.R --sample_set $$*")

sufam/%.tsv : sufam/%.txt
	$$(call RUN,-s 4G -m 6G,"$(RSCRIPT) modules/variant_callers/updatesamples.R --sample_set $$*")
	
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples,$(set))))
		
endif

summary/sufam_summary.xlsx : $(wildcard $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv))
	$(call RUN,-s 12G -m 16G,"export R_LIBS='~/share/usr/anaconda-envs/jrflab-modules-0.1.4/lib/R/library:~/share/usr/lib64/R/library' && \
							  $(RSCRIPT) modules/summary/sufamsummary.R --sample_sets '$(SAMPLE_SETS)'")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)