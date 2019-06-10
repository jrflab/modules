include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/qdnaseq_extract.$(NOW)
PHONY += qdnaseq qdnaseq/readcounts qdnaseq/isobars qdnaseq/variance qdnaseq/log2ratio qdnaseq/bed

qdnaseq_extract : $(foreach sample,$(SAMPLES),qdnaseq/readcounts/$(sample).pdf qdnaseq/isobars/$(sample).pdf qdnaseq/variance/$(sample).pdf qdnaseq/log2ratio/$(sample).pdf qdnaseq/bed/$(sample).bed)

DEFAULT_ENV = $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.6
QDNASEQ_ENV = $(HOME)/share/usr/anaconda-envs/qdnaseq
QDNASEQ_BINSIZE = 5

define qdnaseq-log2ratio
qdnaseq/readcounts/%.pdf qdnaseq/isobars/%.pdf qdnaseq/variance/%.pdf qdnaseq/log2ratio/%.pdf qdnaseq/bed/%.bed : bam/%.bam
	$$(call RUN,-c -n 16 -s 2G -m 3G -w 7200 -v $$(DEFAULT_ENV),"source activate $$(QDNASEQ_ENV) && \
																 $$(RSCRIPT) modules/test/copy_number/qdnaseqextract.R --sample $$(*) --binsize $(QDNASEQ_BINSIZE)")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call qdnaseq-log2ratio,$(sample))))
	
.PHONY: $(PHONY)
