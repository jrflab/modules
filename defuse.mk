# This module is used for running defuse
# input: $(SAMPLES) 
# Options: PHRED64 = true/false
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
#
include ~/share/modules/Makefile.inc

SAMPLE_FILE ?= samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

DEFUSE_CONFIG_FILE = $(HOME)/share/usr/defuse-0.6.1/scripts/config.txt

LOGDIR = log

# Runs defuse locally on the same node
LOCAL ?= FALSE

# Only applies if LOCAL is set to TRUE
NUM_CORES ?= 2

ifeq ($(LOCAL),true)
	DEFUSE_OPTS = -p $(NUM_CORES)
else
	DEFUSE_OPTS = -s sge
endif

.PHONY : defuse

defuse : $(foreach sample,$(SAMPLES),defuse/$(sample).defuse.timestamp)

defuse/%.defuse.timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call INIT_MEM,1G,2G) \
	$(DEFUSE) -c $(DEFUSE_CONFIG_FILE) -1 $(word 1,$^) -2 $(word 2,$^) -o $(@D)/$* $(DEFUSE_OPTS) &> $(LOG) && touch $@
