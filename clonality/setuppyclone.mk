include modules/Makefile.inc

LOGDIR ?= log/setup_pyclone.$(NOW)
PHONY += pyclone

setup_pyclone : $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/config.yaml)

define make-input-pyclone
pyclone/%/config.yaml : sufam/%.tsv
	$$(call RUN, -s 4G -m 6G,"mkdir -p pyclone/$$(*) && \
							  $(RSCRIPT) modules/clonality/tsvforpyclone.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) && \
							  $(RSCRIPT) modules/clonality/pycloneconfig.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES)")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call make-input-pyclone,$(set))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
