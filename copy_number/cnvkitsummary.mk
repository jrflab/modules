include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_summary.$(NOW)
PHONY += cnvkit cnvkit/summary

cnvkit_summary : cnvkit/summary/bygene.txt cnvkit/summary/bygene.RData

facets/summary/bygene.txt facets/summary/bygene.RData : $(foreach pair,$(TUMOR_SAMPLES),cnvkit/called/$(sample).RData)
	$(call RUN,-c -s 8G -m 30G,"mkdir -p cnvkit/summary && \
							 	$(RSCRIPT) modules/copy_number/cnvkitsummary.R --sample_names '$(TUMOR_SAMPLES)'")
												 
.PHONY: $(PHONY)
