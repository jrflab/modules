# Run fusioncatcher
##### DEFAULTS ######

REF ?= hg19
LOGDIR = log/fusioncatcher.$(NOW)

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

FUSIONCATCHER = $(HOME)/share/usr/fusioncatcher/bin/fusioncatcher
FUSIONCATCHER_OPTS = -d $(HOME)/share/usr/fusioncatcher/data/current --extract-buffer-size=10000000000

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

all : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample).fusioncatcher_timestamp)

fusioncatcher/%.fusioncatcher_timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,8,1G,3G,"$(FUSIONCATCHER) $(FUSIONCATCHER_OPTS) -p 8 -o $(@D)/$* -i $<$(,)$(<<) && touch $@")
