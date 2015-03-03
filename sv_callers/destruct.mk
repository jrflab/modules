# Run destruct
#
#
# Author: Fong Chun Chan <fongchunchan@gmail.com>
#

SHELL := /bin/bash

include modules/Makefile.inc

LOGDIR = destruct/log

DESTRUCT_CONFIG_FILE = $(HOME)/share/usr/destruct/destruct/config.txt
#ANALYZE_DNA_BAM = $(HOME)/share/usr/nfuse-0.1.2/scripts/analyze_dna_bam.pl -c $(HOME)/usr/nfuse-0.1.2/scripts/config.txt
#DESTRUCT = /scratch/sohrab_temp/amcpherson_tmp/forray/install/bin/python2.7 /scratch/sohrab_temp/amcpherson_tmp/forray/install/bin/destruct.py /scratch/sohrab_temp/amcpherson_tmp/forray/genesis_config.ini
#DESTRUCT = /scratch/sohrab_temp/amcpherson_tmp/forray/install/bin/python2.7 $(HOME)/share/usr/destruct/destruct/destruct.py /scratch/sohrab_temp/amcpherson_tmp/config.ini
DESTRUCT = $(PYTHON) $(HOME)/share/usr/destruct/destruct/destruct.py $(DESTRUCT_CONFIG_FILE)

VPATH = bam

SAMPLE_FILE = samples.txt

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all

all : destruct/all.timestamp

####
## Build the destruct file list
####
destruct/all_file_list.txt : $(foreach sample,$(SAMPLES),$(sample).bam)
	mkdir -p $(@D); rm -f $@; for bam in $^; do \
		sample=`echo "$$bam" | sed 's/.*\///; s/\..*//;'`; \
		echo -e "$$sample\t$$bam" >> $@; \
	done
	
destruct/%.timestamp : destruct/%_file_list.txt
	mkdir -p $(@D)/$*.tmp $(@D)/breakpoints $(@D)/breakreads $(LOGDIR); $(DESTRUCT) $< $(@D)/$*.tmp $(@D)/breakpoints/$*.breakpoints $(@D)/breakreads/$*.breakreads qsub -p 100 &> $(LOGDIR)/$*.log && touch $@
