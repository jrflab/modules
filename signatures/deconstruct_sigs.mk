include modules/Makefile.inc

LOGDIR = log/deconstruct_sigs.$(NOW)
PHONY += deconstructsigs deconstructsigs/signatures deconstructsigs/plots

deconstructsigs : $(foreach sample,$(TUMOR_SAMPLES),deconstructsigs/signatures/$(sample).RData) $(foreach sample,$(TUMOR_SAMPLES),deconstructsigs/plots/trint_context/$(sample).pdf)

define extract-signatures
deconstructsigs/signatures/%.RData : summary/tsv/mutation_summary.tsv
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"$(RSCRIPT) modules/signatures/extract_signatures.R --sample_name $$(*)")
	
deconstructsigs/plots/trint_context/%.pdf : deconstructsigs/signatures/%.RData
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"mkdir -p  deconstructsigs/plots/trint_context && \
													   mkdir -p  deconstructsigs/plots/signature_exposures && \
													   $(RSCRIPT) modules/signatures/plot_signatures.R --sample_name $$(*)")

endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call extract-signatures,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
