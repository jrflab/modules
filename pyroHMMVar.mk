# run pyrohmmvar: realignment-based variant calling method for 454 and ion torrent

include ~/share/modules/Makefile.inc

LOGDIR = log/pyrohmm.$(NOW)

PYROHMMVAR = $(HOME)/share/usr/bin/pyrohmmvar
PYROHMMVAR_MODEL = $(HOME)/share/reference/pyrohmm_parameter_config

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: all

all : $(foreach sample,$(SAMPLES),pyrohmm/tables/$(sample).pyrohmm.txt)

pyrohmm/tables/%.pyrohmm.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,6G,10G,"$(PYROHMMVAR) -b $< -f $(REF_FASTA) -m $(PYROHMMVAR_MODEL) > $@")


