include modules/Makefile.inc

LOGDIR ?= log/qdnaseq.$(NOW)
PHONY += qdnaseq

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)