include modules/Makefile.inc

LOGDIR ?= log/ms_pyclone.$(NOW)
PHONY += pyclone sufam summary pyclone

PYCLONE_WORKFLOW += sufam_multisample
PYCLONE_WORKFLOW += setup_pyclone
PYCLONE_WORKFLOW += run_pyclone
PYCLONE_WORKFLOW += plot_pyclone

pyclone_workflow : $(PYCLONE_WORKFLOW)

include modules/variant_callers/sufammultisample.mk
include modules/clonality/setuppyclone.mk
include modules/clonality/runpyclone.mk
include modules/clonality/plotpyclone.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
