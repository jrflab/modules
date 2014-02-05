# This module is used for running defuse
# input: $(SAMPLES) 
# Options: PHRED64 = true/false
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
#
include ~/share/modules/Makefile.inc

DEFUSE_CONFIG_FILE = $(HOME)/share/usr/defuse-0.6.1/scripts/config.txt
#DEFUSE_CONFIG_FILE = /opt/common/defuse/defuse-0.6.1/scripts/config.txt
DEFUSE_FILTER = $(HOME)/share/scripts/filterDefuse.pl
DEFUSE_NORMAL_FILTER = $(HOME)/share/scripts/normalFilterDefuse.pl

LOGDIR = log/defuse.$(NOW)

# Runs defuse locally on the same node
LOCAL ?= FALSE

# Only applies if LOCAL is set to TRUE
NUM_CORES ?= 2

ifeq ($(LOCAL),true)
	DEFUSE_OPTS = -p $(NUM_CORES)
else
	DEFUSE_OPTS = -s sge -p 10
endif

.PHONY : all tables


#all : $(foreach sample,$(SAMPLES),defuse/$(sample).defuse_timestamp)
ifdef NORMAL_DEFUSE_RESULTS
all : defuse/tables/all.defuse_results.nft.txt tables
else
all : defuse/tables/all.defuse_results.txt tables
endif

tables : $(foreach sample,$(SAMPLES),defuse/tables/$(sample).defuse_results.txt)

defuse/%.defuse_timestamp : fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5
	$(INIT) $(CHECK_MD5) $(DEFUSE) -c $(DEFUSE_CONFIG_FILE) -1 $(word 1,$(^M)) -2 $(word 2,$(^M)) -o $(@D)/$* $(DEFUSE_OPTS) &> $(LOG) && touch $@

defuse/tables/%.defuse_results.txt : defuse/%.defuse_timestamp
	$(INIT) $(PERL) $(DEFUSE_FILTER) defuse/$*/results.filtered.tsv > $@ 2> $(LOG) && rm -r defuse/$*

defuse/tables/all.defuse_results.txt : $(foreach sample,$(SAMPLES),defuse/tables/$(sample).defuse_results.txt)
	$(INIT) head -1 $< > $@ && for x in $^; do sed '1d' $$x >> $@; done

defuse/tables/%.defuse_results.nft.txt : defuse/tables/%.defuse_results.txt $(NORMAL_DEFUSE_RESULTS)
	$(INIT) $(PERL) $(DEFUSE_NORMAL_FILTER) -w 1000 $(NORMAL_DEFUSE_RESULTS) $< > $@

include ~/share/modules/fastq.mk
