include modules/Makefile.inc

LOGDIR = log/sv_signature.$(NOW)

svsignature : $(foreach sample,$(TUMOR_SAMPLES),sv_signature/$(sample)/)

define extract-signatures
sv_signature/$1/ : 
	$$(call RUN,-s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"set -o pipefial && \
							   $(RSCRIPT) modules/signatures/extract_signatures.R \
							   --sample_name $$()")
	
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call extract-signatures,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(DECONSTRUCTSIGS_ENV)/bin/R --version > version/deconstruct_sigs.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: deconstructsigs