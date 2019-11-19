include modules/Makefile.inc

LOGDIR ?= log/delmh_summary.$(NOW)
PHONY += delmh_summary

delmh_summary : summary/tsv/delmh_summary.tsv

summary/tsv/delmh_summary.tsv : summary/tsv/mutation_summary.tsv
	$(call RUN,-n 1 -s 8G -m 8G,"set -o pipefail && \
								 $(RSCRIPT) modules/summary/delmh_summary.R --input_file $(<)")
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
