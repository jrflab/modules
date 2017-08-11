# Run fusionfinder
##### DEFAULTS ######

LOGDIR = log/fusionfinder.$(NOW)

##### MAKE INCLUDES #####
include modules/Makefile.inc

FUSIONFINDER = 
FUSIONCATCHER_OPTS = --phred33 --cref  --ncref

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

all : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample).fusioncatcher_timestamp)

fusioncatcher/%.fusioncatcher_timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-n 8 -s 1G -m 4G,"$(FUSIONCATCHER) $(FUSIONCATCHER_OPTS) -p 8 -o $(@D)/$* -i $<$(,)$(<<) && touch $@")
