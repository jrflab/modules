include modules/Makefile.inc

LOGDIR = log/deconstruct_sigs.$(NOW)
PHONY += deconstructsigs deconstructsigs/signatures deconstructsigs/plots

deconstructsigs : $(foreach sample,$(TUMOR_SAMPLES),deconstructsigs/signatures/$(sample).RData) $(foreach sample,$(TUMOR_SAMPLES),deconstructsigs/plots/trint_context/$(sample).pdf)

SUFAM ?= false

ifeq ($(SUFAM),true)

define extract-signatures
deconstructsigs/signaures/%.RData : summary/tsv/mutation_summary.tsv
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"$(RSCRIPT) modules/mut_sigs/extract_signatures.R --sample_name $$(*)")
	
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call extract-signatures,$(sample))))
		
else 

define extract-signatures
deconstructsigs/signatures/%.RData : summary/tsv/mutation_summary.tsv
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"$(RSCRIPT) modules/mut_sigs/extract_signatures.R --sample_name $$(*)")
	
deconstructsigs/plots/trint_context/%.pdf : deconstructsigs/signatures/%.RData
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"mkdir -p  deconstructsigs/plots/trint_context && \
													   mkdir -p  deconstructsigs/plots/signature_exposures && \
													   $(RSCRIPT) modules/mut_sigs/plot_signatures.R --sample_name $$(*)")

endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call extract-signatures,$(sample))))
		
endif


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
