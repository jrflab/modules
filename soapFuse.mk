# soapfuse
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR = log/soapfuse.$(NOW)

SOAPFUSE = $(HOME)/usr/SOAPfuse-v1.26/SOAPfuse-RUN.pl 
SOAPFUSE_CONFIG = $(HOME)/share/usr/SOAPfuse-v1.26/config/config.txt
PREPARE_SOAPFUSE = $(HOME)/share/scripts/prepareSoapFuse.pl

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all sample_tables soapfuse

all : soapfuse/tables/all.sfuse.txt soapfuse/tables/all.sfuse.isoforms.txt
sample_tables : $(foreach sample,$(SAMPLES),soapfuse/tables/$(sample).sfuse.txt soapfuse/tables/$(sample).sfuse.isoforms.txt)
soapfuse : $(foreach sample,$(SAMPLES),soapfuse/$(sample).timestamp)

soapfuse/%.timestamp : soapfuse/sample_lists/%.txt
	$(call LSCRIPT_NAMED_PARALLEL_MEM,$*_soapfuse,4,2G,2.5G,"$(SOAPFUSE) -c $(SOAPFUSE_CONFIG) -fd $(@D) -l $< -o $(@D)/$* && touch $@")

soapfuse/sample_lists/%.txt : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(INIT) SAMPLEDIR=soapfuse/$*/$*; \
	mkdir -p $$SAMPLEDIR; ln -f -t $$SAMPLEDIR $^; \
	rename .1.fastq.gz _1.fastq.gz $$SAMPLEDIR/*.gz; \
	rename .2.fastq.gz _2.fastq.gz $$SAMPLEDIR/*.gz; \
	READLEN=`zcat $< | head -2 | sed 1d | wc -c`; \
	echo -e "$*\t$*\t$*\t$$READLEN" > $@

soapfuse/tables/%.sfuse.txt : soapfuse/%.timestamp
	$(INIT) ln -t $(@D) soapfuse/$*/final_fusion_genes/$*/$*.final.Fusion.specific.for.genes

soapfuse/tables/%.sfuse.isoforms.txt : soapfuse/%.timestamp
	$(INIT) ln -t $(@D) soapfuse/$*/final_fusion_genes/$*/$*.final.Fusion.specific.for.trans

soapfuse/tables/all.sfuse.%.txt : $(foreach sample,$(SAMPLES),soapfuse/tables/$(sample).sfuse.%.txt)
	$(INIT) head -1 $< | sed 's/^/Sample\t/;' > $@ && for i in $^; do sed "1d; s/^/$$(basename $${i%%.sfuse.$*.txt})\t/" $$i >> $@; done
