include modules/Makefile.inc

LOGDIR ?= log/qdnaseq.$(NOW)
PHONY += qdnaseq qdnaseq/readcounts qdnaseq/isobars qdnaseq/variance qdnaseq/log2ratio qdnaseq/bed

qdnaseq : $(foreach sample,$(SAMPLES),qdnaseq/readcounts/$(sample).pdf qdnaseq/isobars/$(sample).pdf qdnaseq/variance/$(sample).pdf qdnaseq/log2ratio/$(sample).pdf qdnaseq/bed/$(sample).bed)

define qdnaseq-log2ratio
qdnaseq/readcounts/%.pdf qdnaseq/isobars/%.pdf qdnaseq/variance/%.pdf qdnaseq/log2ratio/%.pdf qdnaseq/bed/%.bed : bam/%.bam
	$$(call RUN,-c -s XXG -m XXG -w 7200,"source $(HOME)/share/usr/opt/miniconda/bin/activate /home/brownd7/share/usr/anaconda-envs/jrflab-modules-0.1.5 && \
										  export CPATH=$(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.5/include:$(CPATH) && \
										  export LIBRARY_PATH='$(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.5/lib:$(LIBRARY_PATH)' && \
										  export LD_LIBRARY_PATH='$(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.5/lib:$(LD_LIBRARY_PATH)' && \
										  export R_LIBS='$(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.5/usr/R/library:$(R_LIBS)' && \
										  export R_LIBS='$(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.5/lib/R/library:$(R_LIBS)' && \
										  $(RSCRIPT) modules/test/copy_number/qdnaseq.R --sample $$(*)")
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)



               
                
                
                
                