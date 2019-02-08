include modules/Makefile.inc
include modules/config.inc

LOGDIR = log/copynumber_summary.$(NOW)
PHONY += genome_stats summary summary/tsv

CN_SUMMARY_WORKFLOW += genome_altered
CN_SUMMARY_WORKFLOW += lst_score
CN_SUMMARY_WORKFLOW += ntai_score
CN_SUMMARY_WORKFLOW += myriad_score
CN_SUMMARY_WORKFLOW += genome_summary

cn_summary_workflow : $(CN_SUMMARY_WORKFLOW)

include modules/copy_number/genomealtered.mk
include modules/copy_number/lstscore.mk
include modules/copy_number/ntaiscore.mk
include modules/copy_number/myriadhrdscore.mk
include modules/summary/genomesummary.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
