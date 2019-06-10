include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/qdnaseq.$(NOW)
PHONY += qdnaseq qdnaseq/copynumber qdnaseq/copynumber/log2ratio qdnaseq/copynumber/segmented qdnaseq/copynumber/pcf

QDNA_SEQ_WORKFLOW += qdnaseq_extract
QDNA_SEQ_WORKFLOW += qdnaseq_copynumber

qdna_seq_workflow : $(QDNA_SEQ_WORKFLOW)

include modules/test/copy_number/qdnaseqextract.mk
include modules/test/copy_number/qdnaseqcopynumber.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
