include modules/Makefile.inc

LOGDIR ?= log/genome_summary.$(NOW)
PHONY += genome_summary

LST_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).lst))
GENOME_ALTERED ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).fga))
NTAI_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).ntai))
MYRIAD_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).mrs))

genome_summary : genome_stats/lst_score.tsv genome_stats/genome_altered.tsv genome_stats/ntai_score.tsv genome_stats/myriad_score.tsv

genome_stats/lst_score.tsv genome_stats/genome_altered.tsv genome_stats/ntai_score.tsv genome_stats/myriad_score.tsv genome_stats/genome_summary.tsv : 
	$(call RUN,-n 1 -s 4G -m 4G,"cat $(LST_SCORE) > genome_stats/lst_score.tsv && \
				     			 cat $(GENOME_ALTERED) > genome_stats/genome_altered.tsv && \
				     			 cat $(NTAI_SCORE) > genome_stats/ntai_score.tsv && \
				     			 cat $(MYRIAD_SCORE) > genome_stats/myriad_score.tsv && \
				     			 $(RSCRIPT) modules/summary/genomesummary.R")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
