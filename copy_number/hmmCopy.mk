include ~/share/modules/Makefile.inc

LOGDIR = log/hmmCopy.$(NOW)

WINDOW_SIZE ?= 1000

READ_COUNTER = $(HOME)/share/usr/bin/readCounter

MAP_BW = $(HOME)/share/references/genomes/wgEncodeCrgMapabilityAlign100mer.bigWig

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all

all : $(foreach sample,$(SAMPLES),hmmCopy/w$(WINDOW_SIZE)_wig/$(sample).wig)

hmmCopy/%.w$(WINDOW_SIZE).wig : bam/%.bam
	$(call LSCRIPT_MEM,4G,7G,"$(READ_COUNTER) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $< > $@")
