include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/qdnaseq_copynumber.$(NOW)
PHONY += qdnaseq qdnaseq/copynumber qdnaseq/copynumber/log2ratio qdnaseq/copynumber/segmented qdnaseq/copynumber/pcf

qdnaseq : $(foreach sample,$(SAMPLES),qdnaseq/copynumber/log2ratio/$(sample).pdf qdnaseq/copynumber/segmented/$(sample).RData qdnaseq/copynumber/pcf/$(sample).pdf)

define qdnaseq-plot-log2ratio
qdnaseq/copynumber/log2ratio/%.pdf : qdnaseq/bed/%.bed
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/qdnaseqplot.R --sample $$(*) --type 'raw'")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call qdnaseq-plot-log2ratio,$(sample))))
		
define qdnaseq-segment-log2ratio
qdnaseq/copynumber/segmented/%.RData qdnaseq/copynumber/pcf/%.pdf : qdnaseq/bed/%.bed
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/qdnaseqsegment.R --sample $$(*)")
	
qdnaseq/copynumber/pcf/%.pdf : qdnaseq/copynumber/segmented/%.RData
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/qdnaseqplot.R --sample $$(*) --type 'segmented'")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call qdnaseq-segment-log2ratio,$(sample))))
		

.PHONY: $(PHONY)
