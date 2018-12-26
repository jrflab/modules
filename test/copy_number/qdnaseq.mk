include modules/Makefile.inc

LOGDIR ?= log/qdnaseq.$(NOW)
PHONY += qdnaseq qdnaseq/readcounts qdnaseq/isobars qdnaseq/variance qdnaseq/log2ratio qdnaseq/bed

qdnaseq : $(foreach sample,$(SAMPLES),qdnaseq/readcounts/$(sample).pdf qdnaseq/isobars/$(sample).pdf qdnaseq/variance/$(sample).pdf qdnaseq/log2ratio/$(sample).pdf qdnaseq/bed/$(sample).bed)

DEFAULT_ENV = $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.6
QDNASEQ_ENV = $(HOME)/share/usr/anaconda-envs/qdnaseq

define qdnaseq-log2ratio
qdnaseq/readcounts/%.pdf qdnaseq/isobars/%.pdf qdnaseq/variance/%.pdf qdnaseq/log2ratio/%.pdf qdnaseq/bed/%.bed : bam/%.bam
	$$(call RUN,-c -n 16 -s 1G -m 3G -w 7200 --default_env $$(DEFAULT_ENV) -v $$(DEFAULT_ENV),"$(DEFAULT_ENV)/bin/activate $(QDNASEQ_ENV) && \
																							   export RLIBS=/lila/data/reis-filho/usr/anaconda-envs/qdnaseq/lib/R/library && \
																							   $(RSCRIPT) modules/test/copy_number/qdnaseq.R --sample $$(*) --binsize $(QDNASEQ_BINSIZE)")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call qdnaseq-log2ratio,$(sample))))
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
