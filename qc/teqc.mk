# Run TEQC R library on bams
# vim: set ft=make :
include modules/Makefile.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))
LOGDIR ?= teqc/log
TARGETS_BED_FILE ?= intervals.bed

.PHONY: teqc
.DELETE_ON_ERROR:
.SECONDARY:

teqc : teqc_report/index.html

teqc/%.Rdata : bam/%.bam bam/%.bam.bai
	$(call INIT_MEM,12G,14G) $(RSCRIPT) $(HOME)/share/scripts/TEQC.R --ref=$(REF) --outFile $@ $< $(TARGETS_BED_FILE) &> $(LOGDIR)/$(@F).log

teqc_report/index.html : $(foreach sample,$(SAMPLES),teqc/$(sample).Rdata)
	$(call INIT_MEM,12G,14G) $(MKDIR) teqc_report; $(RSCRIPT) $(HOME)/share/scripts/TEQCreport.R --outDir=$(@D) $^

