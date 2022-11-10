include modules/Makefile.inc

LOGDIR = log/deconstruct_sigs.$(NOW)

deconstructsigs : $(foreach sample,$(TUMOR_SAMPLES),deconstructsigs/signatures/$(sample).RData) \
		  $(foreach sample,$(TUMOR_SAMPLES),deconstructsigs/plots/context/$(sample).pdf)

define extract-signatures
deconstructsigs/signatures/%.RData : summary/tsv/mutation_summary.tsv
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"set -o pipefial && \
							   $(RSCRIPT) modules/signatures/extract_signatures.R \
							   --sample_name $$()")
	
deconstructsigs/plots/context/%.pdf : deconstructsigs/signatures/%.RData
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"set -o pipefail && \
							   $(RSCRIPT) modules/signatures/plot_signatures.R \
							   --sample_name $$(*)")

endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call extract-signatures,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(DECONSTRUCTSIGS_ENV)/bin/R --version > version/deconstruct_sigs.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: deconstructsigs