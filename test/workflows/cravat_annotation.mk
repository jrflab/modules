include modules/Makefile.inc
include modules/config.inc

LOGDIR = log/cravat_annotation.$(NOW)
PHONY += cravat summary summary/tsv

ANNOTATION_WORKFLOW += cravat_annotate
ANNOTATION_WORKFLOW += cravat_summary

cravat_annotation_workflow : $(ANNOTATION_WORKFLOW)

include modules/vcf_tools/cravat_annotation.mk
include modules/summary/cravat_summary.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


