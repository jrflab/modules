include modules/Makefile.inc

LOGDIR ?= log/msk_access.$(NOW)
PHONY += bam fgbio

MSK_ACCESS_WORKFLOW += collapsed_umi
MSK_ACCESS_WORKFLOW += mearge_alignments
MSK_ACCESS_WORKFLOW += call_consensus
MSK_ACCESS_WORKFLOW += align_consensus

msk_access_workflow : $(MSK_ACCESS_WORKFLOW)

include modules/test/bam_tools/umicollapsing.mk
include modules/test/bam_tools/mergealignments.mk
include modules/test/bam_tools/callconsensus.mk
include modules/test/bam_tools/alignconsensus.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
