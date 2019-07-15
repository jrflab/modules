include modules/Makefile.inc

LOGDIR ?= log/marianas_access.$(NOW)
PHONY += bam marianas cnvaccess

MSK_ACCESS_WORKFLOW += clip_umi
MSK_ACCESS_WORKFLOW += align_fastq

msk_access_workflow : $(MSK_ACCESS_WORKFLOW)

include modules/test/fastq_tools/clipumi.mk
include modules/test/bam_tools/alignfastq.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
