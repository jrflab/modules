include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/qdnaseq_copynumber.$(NOW)
PHONY += qdnaseq qdnaseq/copynumber qdnaseq/copynumber/log2ratio

plot : $(foreach sample,$(SAMPLES),qdnaseq/copynumber/log2ratio/$(sample).pdf)

define qdnaseq-plot-log2ratio
qdnaseq/copynumber/log2ratio/%.pdf : qdnaseq/bed/%.bed
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/qdnaseqplot.R --sample $$(*)")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call qdnaseq-plot-log2ratio,$(sample))))
				
.PHONY: $(PHONY)
