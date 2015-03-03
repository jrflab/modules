# Run nfuse on dna bams
# This module is defunct now. Use destruct to call rearrangements for DNA only
#
#
# Author: Raymond Lim <raylim@mm.st>
#

SHELL := /bin/bash

include modules/Makefile.inc

LOGDIR = nfuse/log
ANALYZE_DNA_BAM = $(HOME)/share/usr/nfuse-0.1.2/scripts/analyze_dna_bam.pl -c $(HOME)/usr/nfuse-0.1.2/scripts/config.txt

ANALYZE_DNA_BAM = /scratch/sohrab_temp/amcpherson_tmp/forray/install/bin/python2.7 /scratch/sohrab_temp/amcpherson_tmp/forray/install/bin/destruct.py /scratch/sohrab_temp/amcpherson_tmp/forray/genesis_config.ini

VPATH = bam

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all

all : nfuse/timestamp

nfuse/file_list.txt : $(foreach sample,$(SAMPLES),$(sample).bam}
	mkdir -p $(@D); rm -f $@; for bam in $^; do \
		sample=`echo "$$bam" | sed 's/.*\///; s/\..*//;'`; \
		echo -e "$$sample\t$$bam" >> $@; \
	done
	
nfuse/timestamp : nfuse/file_list.txt
	mkdir -p tmp $(LOGDIR); $(ANALYZE_DNA_BAM) $< tmp $(@D) sge -p 100 &> $(LOGDIR)/nfuse.log && touch $@

