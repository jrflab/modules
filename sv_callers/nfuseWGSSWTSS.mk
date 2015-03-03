# vim: set ft=make :

include modules/Makefile.inc

WGSS_WTSS_PAIR_FILE ?= wgss_wtss_pairs.txt

WGSS_SAMPLES ?= $(shell cut -f 1 $(WGSS_WTSS_PAIR_FILE))
WTSS_SAMPLES ?= $(shell cut -f 2 $(WGSS_WTSS_PAIR_FILE))
SAMPLES ?= $(WGSS_SAMPLES) $(WTSS_SAMPLES)
NSAMPLES ?= $(words $(WGSS_SAMPLES))

$(foreach i,$(shell seq 1 $(NSAMPLES)),$(eval wgss_lookup.$(word $i,$(WGSS_SAMPLES)) := $(word $i,$(WTSS_SAMPLES))))

NFUSE = $(HOME)/share/usr/nfuse-0.2.0/scripts/nfuse.pl -c $(HOME)/share/usr/nfuse-0.2.0/scripts/config.txt -s sge -p 100

#VPATH = ../WTSS/bam ../WGSS/bam

LOGDIR ?= log/nfuse.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach wgss_sample,$(WGSS_SAMPLES),nfuse/$(wgss_sample)_$(wgss_lookup.$(wgss_sample)).timestamp)

#all : $(foreach i,$(shell seq 1 $(NSAMPLES)),nfuse/$(word $i,$(WGSS_SAMPLES))_$(word $i,$(WTSS_SAMPLES)).timestamp)

#$(call nfuse-wgss-wtss,wgss-sample,wtss-sample)
define nfuse-wgss-wtss
nfuse/$1_$2.timestamp : fastq/$1.1.fastq.gz fastq/$1.2.fastq.gz fastq/$2.1.fastq.gz fastq/$2.2.fastq.gz
	$$(INIT) $$(NFUSE) --dnafq1 $$(word 1,$$^) --dnafq2 $$(word 2,$$^) --rnafq1 $$(word 3,$$^) --rnafq2 $$(word 4,$$^) -o nfuse/$1_$2 -n $1_$2 && touch $$@ &> $(LOGDIR)/$1_$2.log
endef
$(foreach i,$(shell seq 1 $(NSAMPLES)),$(eval $(call nfuse-wgss-wtss,$(word $i,$(WGSS_SAMPLES)),$(word $i,$(WTSS_SAMPLES)))))


#include modules/fastq_tools/fastq.mk
