include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc
include ~/share/modules/hg19.inc

LOGDIR = hmmCopy/log

WINDOW_SIZE ?= 1000

READ_COUNTER = $(HOME)/share/usr/bin/readCounter
MAP_COUNTER = $(HOME)/share/usr/bin/mapCounter
GC_COUNTER = $(HOME)/share/usr/bin/gcCounter

MAP_BW = $(HOME)/share/references/genomes/wgEncodeCrgMapabilityAlign100mer.bigWig

SAMPLE_FILE = samples.txt
#SAMPLE_PAIR_FILE = sample_pairs.txt

#TUMOR_SAMPLES := $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
#NORMAL_SAMPLES := $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES = $(shell cat $(SAMPLE_FILE))

VPATH = bam

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all

all : $(foreach sample,$(SAMPLES),hmmCopy/w$(WINDOW_SIZE)_wig/$(sample).wig) hmmCopy/w$(WINDOW_SIZE)_wig/map.wig hmmCopy/w$(WINDOW_SIZE)_wig/gc.wig

hmmCopy/wig/map.wig :
	SGE_RREQ="${SGE_RREQ} ${call MEM_FREE,2G,3G}" $(MKDIR) $(@D) $(LOGDIR); $(MAP_COUNTER) -w $(WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(MAP_BW) > $@

hmmCopy/wig/gc.wig :
	SGE_RREQ="${SGE_RREQ} ${call MEM_FREE,2G,3G}" $(MKDIR) $(@D) $(LOGDIR); $(GC_COUNTER) -w $(WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(REF_FASTA) > $@

hmmCopy/w$(WINDOW_SIZE)_wig/%.wig : %.bam
	SGE_RREQ="${SGE_RREQ} ${call MEM_FREE,2G,3G}" $(MKDIR) $(@D) $(LOGDIR); $(READ_COUNTER) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $< > $@ 2> $(LOGDIR)/$(@F).log
