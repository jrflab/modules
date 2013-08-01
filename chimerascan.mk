# Chimerascan
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/chimscan.$(NOW)
SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

all : $(foreach sample,$(SAMPLES),chimerascan/$(sample).chimscan_timestamp)

CHIMERASCAN = $(PYTHON) /ifs/opt/common/python/python-2.7.3/lib/python2.7/site-packages/chimerascan/chimerascan_run.py
CHIMERASCAN_OPTS = -v --quals illumina

chimerascan/%.chimscan_timestamp : fastq/%.1.fastq fastq/%.2.fastq
	$(call LSCRIPT_PARALLEL_MEM,4,2G,3G, "$(CHIMERASCAN) $(CHIMERASCAN_OPTS) -p 4 $(REF_FASTA) $^ $(@D)/$* && touch $@")

