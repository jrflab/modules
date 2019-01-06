include modules/Makefile.inc

LOGDIR ?= log/cravat_summary.$(NOW)
PHONY += cravat summary summary/tsv

cravat : summary/tsv/cravat_summary.tsv

summary/tsv/cravat_summary.tsv : $(wildcard cravat/$(SAMPLES).txt)
	$(call RUN,-c -s 24G -m 48G -w 7200,"$(RSCRIPT) modules/test/annotations/cravat_summary.R --sample_names '$(SAMPLES)'")
	

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)