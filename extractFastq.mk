# This module extract fastq files from a bam file.  It will use either the Picard (SamToFastq.jar) or bam2fastq programs to extract the fastq.  You can specify which program to use with the EXTRACT_TOOL variable
# input: $(SAMPLES)
# Author: Fong Chun Chan <fongchunchan@gmail.com>

include ~/share/modules/Makefile.inc
include ~/share/modules/hg19.inc

SAM_TO_FASTQ = $(JAVA) -Xmx4G -jar $(JARDIR)/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT

SAMPLE_FILE ?= samplesToExtract.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR ?= log/fastq.$(NOW)

.DELETE_ON_ERROR:

.SECONDARY:

EXTRACT_TOOL ?= PICARD

all : $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz)

ifeq (${EXTRACT_TOOL},PICARD)
fastq/%.1.fastq fastq/%.2.fastq : gsc_bam/%.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,10G,40G)" $(MKDIR) $(LOGDIR);\
	$(SAM_TO_FASTQ) I=$< FASTQ=fastq/$*.1.fastq SECOND_END_FASTQ=fastq/$*.2.fastq &> ${LOGDIR}/$(@F).SamToFastq.log && $(RM) $<
else
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : gsc_bam/%.bam
	$(INIT_MEM,10G,40G) $(BAM2FASTQ) -o fastq/$*#.fastq $< &> $(LOGDIR)/$(@F).bam2fastq.log && mv fastq/$*_1.fastq fastq/$*.1.fastq && mv fastq/$*_2.fastq fastq/$*.2.fastq
endif

fastq/%.fastq.gz : fastq/%.fastq
	$(INIT_MEM,1G,2G) gzip $<
