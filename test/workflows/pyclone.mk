include modules/Makefile.inc

LOGDIR ?= log/cnvkit.$(NOW)
PHONY += sufam summary pyclone

PYCLONE_WORKFLOW += cnvkit_coverage


pyclone_workflow : $(PYCLONE_WORKFLOW)

include modules/copy_number/cnvkitcoverage.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
