# soapfuse
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR = log/soapfuse.$(NOW)

SOAPFUSE = $(HOME)/usr/SOAPfuse-v1.26/SOAPfuse-RUN.pl 
SOAPFUSE_CONFIG = $(HOME)/share/usr/SOAPfuse-v1.26/config/config.txt
PREPARE_SOAPFUSE = $(HOME)/share/scripts/prepareSoapFuse.pl

SOAPFUSE_SAMPLES_FILE = samples.soapfuse.txt

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach sample,$(SAMPLES),soapfuse/$(sample).timestamp)

soapfuse/%.timestamp : soapfuse/sample_lists/%.txt
	$(call LSCRIPT_NAMED_PARALLEL_MEM,$*_soapfuse,4,2G,2.5G,"$(SOAPFUSE) -c $(SOAPFUSE_CONFIG) -fd $(@D) -l $< -o $(@D)/$* && touch $@")

soapfuse/sample_lists/%.txt : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(INIT) SAMPLEDIR=soapfuse/$*/$*; \
	mkdir -p $$SAMPLEDIR; ln -f -t $$SAMPLEDIR $^; \
	rename .1.fastq.gz _1.fastq.gz $$SAMPLEDIR/*.gz; \
	rename .2.fastq.gz _2.fastq.gz $$SAMPLEDIR/*.gz; \
	READLEN=`zcat $< | head -2 | sed 1d | wc -c`; \
	echo -e "$*\t$*\t$*\t$$READLEN" > $@

