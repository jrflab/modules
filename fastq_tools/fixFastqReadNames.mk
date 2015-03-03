# This module is used for fixing read names of paired fastq files
# input: $(SAMPLES)
# Author: Fong Chun Chan <fongchunchan@gmail.com>

include modules/Makefile.inc
include modules/hg19.inc

FIX_FASTQ_READ_NAMES = $(PYTHON) ~/share/scripts/fixFastqReadNames.py

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR = fastq/logs

.DELETE_ON_ERROR:

.SECONDARY:

all : $(foreach sample,${SAMPLES},fastq/$(sample).1.fixed.fastq)

fastq/%.1.fixed.fastq fastq/%.2.fixed.fastq : fastq/%.1.fastq fastq/%.2.fastq
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,2G,5G)" $(MKDIR) $(LOGDIR);\
	$(FIX_FASTQ_READ_NAMES) $(word 1,$^) $(word 2,$^) fastq/$*.1.fixed.fastq fastq/$*.2.fixed.fastq &> ${LOGDIR}/$(@F).log;\
