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
CHIMSCAN_NORMAL_FILTER = $(HOME)/share/scripts/normalFilterChimerascan.pl

RECURRENT_FUSIONS = $(RSCRIPT) $(HOME)/share/scripts/recurrentFusions.R

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

ALL = $(foreach sample,$(SAMPLES),chimscan/$(sample).chimscan_timestamp)
ifdef NORMAL_CHIMSCAN_RESULTS
ALL += $(foreach sample,$(SAMPLES),chimscan/tables/$(sample).chimscan_results.nft.txt)
ALLTABLE = chimscan/alltables/all.chimscan_results.nft.txt
else 
ALLTABLE = chimscan/alltables/all.chimscan_results.txt
endif

ALL += $(ALLTABLE)
ALL += chimscan/recur_tables/recurGenes.txt

all : $(ALL)

CHIMERASCAN = PYTHONPATH=$(CHIMSCAN_PYTHONPATH) $(CHIMSCAN_PYTHON) /home/limr/share/usr/lib/python/chimerascan/chimerascan_run.py
CHIMERASCAN_OPTS = -v --quals illumina

chimscan/%.chimscan_timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,4,6G,12G,"$(CHIMERASCAN) $(CHIMERASCAN_OPTS) -p 4 $(CHIMSCAN_INDEX) $^ $(@D)/$* && touch $@")

#chimerascan/tables/all.chimscan_results.txt : $(foreach sample,$(SAMPLES),chimerascan/$(sample).chimscan_timestamp)
#$(INIT) head -1 $(basename $<)/chimeras.bedpe > $@ && for x in $(addsuffix /chimeras.bedpe,$(basename $^)); do sed '1d' $$x >> $@; done

chimscan/tables/%.chimscan_results.txt : chimscan/%.chimscan_timestamp
	$(INIT) ln -f $(basename $<)/chimeras.bedpe $@ && rm chimscan/$*

%.chimscan_results.nft.txt : %.chimscan_results.txt $(NORMAL_CHIMSCAN_RESULTS)
	$(call LSCRIPT_MEM,2G,4G,"$(PERL) $(CHIMSCAN_NORMAL_FILTER) -w 1000 $(NORMAL_CHIMSCAN_RESULTS) $< > $@")

chimscan/alltables/all.chimscan_results.txt : $(foreach sample,$(SAMPLES),chimscan/tables/$(sample).chimscan_results.txt)
	head -1 $< | sed 's/^/Sample\t/; s/#//' > $@ && for i in $^; do sed "1d; s/^/$$(basename $${i%%.chimscan_results.txt})\t/" $$i >> $@; done

chimscan/recur_tables/recurGenes.txt : $(ALLTABLE)
	$(INIT) $(RECURRENT_FUSIONS) --geneCol1 genes5p --geneCol2 genes3p --sampleCol Sample --outDir $(@D) $< 
