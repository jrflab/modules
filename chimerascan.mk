# Chimerascan
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/chimscan.$(NOW)
SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

CHIMSCAN_PYTHONPATH := /home/limr/share/usr/lib/python:/home/limr/share/usr/lib/python2.7
CHIMSCAN_PYTHON := $(HOME)/share/usr/bin/python
CHIMSCAN_INDEX := $(HOME)/share/reference/chimerascan_index

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

all : $(foreach sample,$(SAMPLES),chimerascan/$(sample).chimscan_timestamp)

CHIMERASCAN = PYTHONPATH=$(CHIMSCAN_PYTHONPATH) $(CHIMSCAN_PYTHON) /home/limr/share/usr/lib/python/chimerascan/chimerascan_run.py
CHIMERASCAN_OPTS = -v --quals illumina

chimerascan/%.chimscan_timestamp : fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5
	$(call LSCRIPT_PARALLEL_MEM,4,6G,10G,"$(CHECK_MD5) $(CHIMERASCAN) $(CHIMERASCAN_OPTS) -p 4 $(CHIMSCAN_INDEX) $(^M) $(@D)/$* && touch $@ && $(RM) -r $(@D)/$*/tmp")

