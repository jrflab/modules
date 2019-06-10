include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/qdnaseq_copynumber.$(NOW)
PHONY += qdnaseq qdnaseq/copynumber qdnaseq/copynumber/log2ratio qdnaseq/copynumber/segmented qdnaseq/copynumber/pcf

qdnaseq_copynumber : $(foreach sample,$(SAMPLES),qdnaseq/copynumber/log2ratio/$(sample).pdf qdnaseq/copynumber/segmented/$(sample).RData qdnaseq/copynumber/pcf/$(sample).pdf)

define qdnaseq-plot-log2ratio
qdnaseq/copynumber/log2ratio/%.pdf : qdnaseq/bed/%.bed
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 10G -m 12G,"$(RSCRIPT) modules/test/copy_number/qdnaseqplot.R --sample $$(*) --type 'raw'")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call qdnaseq-plot-log2ratio,$(sample))))
		
define qdnaseq-segment-log2ratio
qdnaseq/copynumber/segmented/%.RData : qdnaseq/bed/%.bed
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 12G -m 16G,"$(RSCRIPT) modules/test/copy_number/qdnaseqsegment.R --sample $$(*)")
	
qdnaseq/copynumber/pcf/%.pdf : qdnaseq/copynumber/segmented/%.RData
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 12G -m 16G,"$(RSCRIPT) modules/test/copy_number/qdnaseqplot.R --sample $$(*) --type 'bychromosome' --rho '$${qdnaseq_rho.$1}' --psi '$${qdnaseq_psi.$1}' --gamma '$${qdnaseq_gamma.$1}' && \
																	 $(RSCRIPT) modules/test/copy_number/qdnaseqplot.R --sample $$(*) --type 'segmented' --rho '$${qdnaseq_rho.$1}' --psi '$${qdnaseq_psi.$1}' --gamma '$${qdnaseq_gamma.$1}'")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call qdnaseq-segment-log2ratio,$(sample))))
		

.PHONY: $(PHONY)
