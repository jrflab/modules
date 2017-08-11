# This module extract fastq files from a bam file.  It will use either the Picard (SamToFastq.jar) or bam2fastq programs to extract the fastq.  You can specify which program to use with the EXTRACT_TOOL variable
# input: $(SAMPLES)
# Author: Fong Chun Chan <fongchunchan@gmail.com>

include modules/Makefile.inc

LOGDIR ?= log/extract_fastq.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: extract_fastq

VPATH = rawdata unprocessed_bam

extract_fastq : $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz)

fastq/%.1.fastq.gz fastq/%.2.fastq.gz : %.bam
	$(call RUN,-n 4 -s 4G -m 9G,"$(SAMTOOLS2) sort -T $(<D)/$* -O bam -n -@ 4 -m 6G $< | $(SAMTOOLS2) fastq -f 1 -1 >(gzip -c > fastq/$*.1.fastq.gz) -2 >(gzip -c > fastq/$*.2.fastq.gz) -")
