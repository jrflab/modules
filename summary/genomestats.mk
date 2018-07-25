include modules/Makefile.inc

LOGDIR ?= log/genome_stats.$(NOW)
PHONY += genome_stats

LST_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).lst))
FRACTION_GENOME_ALTERED ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).fga))

genome_stats : genome_stats/summary.lst genome_stats/summary.fga genome_stats/summary.tsv

genome_stats/summary.lst genome_stats/summary.fga genome_stats/summary.tsv : 
	$(call RUN,-n 1 -s 4G -m 4G,"cat $(LST_SCORE) > genome_stats/summary.lst && cat $(FRACTION_GENOME_ALTERED) > genome_stats/summary.fga && cut -f 2 genome_stats/summary.lst | paste genome_stats/summary.fga - > genome_stats/summary.tsv")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
