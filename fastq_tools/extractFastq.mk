# This module extract fastq files from a bam file.  It will use either the Picard (SamToFastq.jar) or bam2fastq programs to extract the fastq.  You can specify which program to use with the EXTRACT_TOOL variable
# input: $(SAMPLES)
# Author: Fong Chun Chan <fongchunchan@gmail.com>

include ~/share/modules/Makefile.inc

SAM_TO_FASTQ = $(JAVA) -Xmx4G -jar $(JARDIR)/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT
BAM_TO_FASTQ = $(HOME)/share/usr/bin/bamToFastq

LOGDIR ?= log/extract_fastq.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: extract_fastq

EXTRACT_TOOL ?= bamToFastq

extract_fastq : $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz.md5)

ifeq (${EXTRACT_TOOL},PICARD)
fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5 : bam/%.bam
	$(call LSCRIPT_MEM,10G,20G,"$(SAM_TO_FASTQ) I=$< FASTQ=>(gzip -c > fastq/$*.1.fastq.gz) SECOND_END_FASTQ=>(gzip -c > fastq/$*.2.fastq.gz) && md5sum fastq/$*.1.fastq.gz > fastq/$*.1.fastq.gz.md5 && md5sum fastq/$*.2.fastq.gz > fastq/$*.2.fastq.gz.md5")
else
fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5 : bam/%.bam
	$(call LSCRIPT_MEM,10G,20G,"$(BAM_TO_FASTQ) -bam $< -fq1 >(gzip -c > fastq/$*.1.fastq.gz) -fq2 >(gzip -c > fastq/$*.2.fastq.gz) \
		&& md5sum fastq/$*.1.fastq.gz > fastq/$*.1.fastq.gz.md5 && md5sum fastq/$*.2.fastq.gz > fastq/$*.2.fastq.gz.md5")
endif
