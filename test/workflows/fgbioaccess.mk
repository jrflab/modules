include modules/Makefile.inc

LOGDIR ?= log/fgbio_access.$(NOW)
PHONY += bam fgbio cnvaccess

MSK_ACCESS_WORKFLOW += collapsed_umi
MSK_ACCESS_WORKFLOW += mearge_alignments
MSK_ACCESS_WORKFLOW += call_consensus
MSK_ACCESS_WORKFLOW += align_consensus
MSK_ACCESS_WORKFLOW += cnvaccess_coverage
MSK_ACCESS_WORKFLOW += cnvaccess_reference
MSK_ACCESS_WORKFLOW += cnvaccess_fix
MSK_ACCESS_WORKFLOW += cnvaccess_plot
MSK_ACCESS_WORKFLOW += cnvaccess_segment

msk_access_workflow : $(MSK_ACCESS_WORKFLOW)

include modules/test/bam_tools/umicollapsing.mk
include modules/test/bam_tools/mergealignments.mk
include modules/test/bam_tools/callconsensus.mk
include modules/test/bam_tools/alignconsensus.mk
include modules/test/copy_number/cnvaccesscoverage.mk
include modules/test/copy_number/cnvaccessreference.mk
include modules/test/copy_number/cnvaccessfix.mk
include modules/test/copy_number/cnvaccessplot.mk
include modules/test/copy_number/cnvaccesssegment.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
