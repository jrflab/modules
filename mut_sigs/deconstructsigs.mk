include modules/Makefile.inc

LOGDIR = log/deconstruct_sigs.$(NOW)
PHONY += deconstructsigs deconstructsigs/signatures

deconstructsigs : $(foreach sample,$(TUMOR_SAMPLES),deconstructsigs/signatures/$(sample).RData)

SUFAM ?= false

ifeq ($(SUFAM),true)

define extract-signatures
deconstructsigs/signaures/%.RData : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"$(RSCRIPT) modules/mut_sigs/ --sample_name $$(*)")
	
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call extract-signatures,$(sample))))
		
else 

define extract-signatures
deconstructsigs/signaures/%.RData : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G -v $(DECONSTRUCTSIGS_ENV),"$(RSCRIPT) modules/mut_sigs/ --sample_name $$(*)")

endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call extract-signatures,$(sample))))
		
endif


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
