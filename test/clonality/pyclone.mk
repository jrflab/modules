include modules/Makefile.inc

LOGDIR ?= log/setup_pyclone.$(NOW)
PHONY += pyclone

pyclone : $(foreach pair,$(SAMPLE_PAIRS),pyclone/$(pair)/config.yaml)

define make-input-pyclone
pyclone/%/config.yaml : summary/tsv/mutation_summary.tsv
	$$(call RUN, -s 16G -m 24G,"mkdir -p pyclone/$$(*) && \
							    $(RSCRIPT) modules/test/clonality/tsvtopyclone.R --sample_name $$(*)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call make-input-pyclone,$(pair))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
