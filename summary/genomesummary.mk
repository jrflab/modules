include modules/Makefile.inc
#include modules/copy_number/genomealtered.mk
include modules/copy_number/lstscore.mk
include modules/copy_number/ntaiscore.mk
#include modules/copy_number/myriadhrdscore.mk

LOGDIR ?= log/genome_summary.$(NOW)

#GENOME_ALTERED = $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).fga))
#LST_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).lst))
#NTAI_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).ntai))
#MYRIAD_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).mrs))

#genome_summary : genome_stats/genome_altered.tsv \
#		 genome_stats/lst_score.tsv \
#		 genome_stats/ntai_score.tsv \
#		 genome_stats/myriad_score.tsv \
#		 summary/tsv/genome_summary.tsv \
#		 summary/genome_summary.xlsx
		 
#genome_summary += genome_altered
#genome_summary += lst_score
#genome_summary += ntai_score
#genome_summary += myriad_score

#genome_stats/genome_altered.tsv : $(GENOME_ALTERED)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $$(GENOME_ALTERED) > $$(@)")
#							 
#genome_stats/lst_score.tsv : $(LST_SCORE)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $(LST_SCORE) > $$(@)")
#				     
#genome_stats/ntai_score.tsv : $(NTAI_SCORE)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $(NTAI_SCORE) > $$(@)")
#
#genome_stats/myriad_score.tsv : $(MYRIAD_SCORE)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $(MYRIAD_SCORE) > $$(@)")
#
#summary/tsv/genome_summary.tsv : genome_stats/genome_altered.tsv genome_stats/lst_score.tsv genome_stats/ntai_score.tsv genome_stats/myriad_score.tsv
#	$(call RUN,-n 1 -s 6G -m 8G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     $(RSCRIPT) modules/summary/genomesummary.R")
#
#summary/genome_summary.xlsx : summary/tsv/genome_summary.tsv
#	$(call RUN,-n 1 -s 4G -m 4G,"python modules/summary/genome_summary_excel.py")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: genome_sumary
