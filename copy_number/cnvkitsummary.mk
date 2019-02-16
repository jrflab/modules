include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_summary.$(NOW)
PHONY += cnvkit cnvkit/summary

cnvkit_summary : cnvkit/summary/bygene.txt

cnvkit/summary/bygene.txt : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/called/$(sample).RData)
	$(call RUN,-c -s 24G -m 48G,"mkdir -p cnvkit/summary && \
							 	 $(RSCRIPT) modules/copy_number/cnvkitsummary.R --sample_names '$(TUMOR_SAMPLES)'")
												 
.PHONY: $(PHONY)
