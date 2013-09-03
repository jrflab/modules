# Run pindel
# Detect indels
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/mutect.$(NOW)
SAMPLE_PAIR_FILE ?= sample_pairs.txt
SAMPLE_FILE ?= samples.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

GET_INSERT_SIZE = $(HOME)/share/usr/bin/getInsertSize.py

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all mutect_vcfs mutect_tables ext_output

all : mutect_vcfs mutect_tables ext_output


