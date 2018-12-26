include modules/Makefile.inc

LOGDIR ?= log/qdnaseq.$(NOW)
PHONY += qdnaseq qdnaseq/readcounts qdnaseq/isobars qdnaseq/variance qdnaseq/log2ratio qdnaseq/bed

qdnaseq : $(foreach sample,$(SAMPLES),qdnaseq/readcounts/$(sample).pdf qdnaseq/isobars/$(sample).pdf qdnaseq/variance/$(sample).pdf qdnaseq/log2ratio/$(sample).pdf qdnaseq/bed/$(sample).bed)

define qdnaseq-log2ratio
qdnaseq/readcounts/%.pdf qdnaseq/isobars/%.pdf qdnaseq/variance/%.pdf qdnaseq/log2ratio/%.pdf qdnaseq/bed/%.bed : bam/%.bam
	$$(call RUN,-c -s XXG -m XXG -w 7200 -v ~/share/usr/anaconda-envs/qdnaseq-1.16.0/,"$(RSCRIPT) modules/test/copy_number/qdnaseq.R --sample $$(*)")
	

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


