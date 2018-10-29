include modules/Makefile.inc

LOGDIR ?= log/extract_unmapped_pairs.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: extract_unmapped_pairs

VPATH = bam
JAVA = $(HOME)/share/usr/jdk1.8.0_121/bin/java
PICARD = /lila/data/riazlab/lib/src/picard.jar

extract_unmapped_pairs : $(foreach sample,$(SAMPLES),extracted_reads/unmapped_pairs/$(sample)_1.fastq)

define extract-unmapped-pairs
extracted_reads/unmapped_pairs/%.bam : extracted_reads/unmapped_pairs/%.txt
	$$(call RUN,-c -n 4 -s 4G -m 9G,"$(JAVA) -jar $(PICARD) FilterSamReads I=unmapped_reads/$$*.bam O=extracted_reads/unmapped_pairs/$$*.bam \
		READ_LIST_FILE=extracted_reads/unmapped_pairs/$$*.txt FILTER=includeReadList")

extracted_reads/unmapped_pairs/%.txt: unmapped_reads/%.bam
	$$(call RUN,-c -n 1 -s 4G -m 9G,"$(SAMTOOLS2) view $$< | cut -f1 | sort | uniq > extracted_reads/unmapped_pairs/$$*.txt")

extracted_reads/unmapped_pairs/%_1.fastq extracted_reads/unmapped_pairs/%_2.fastq : extracted_reads/unmapped_pairs/%.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"bamToFastq -i $$< -fq extracted_reads/unmapped_pairs/$$*_1.fastq -fq2 extracted_reads/unmapped_pairs/$$*_2.fastq")

endef
$(foreach pair,$(SAMPLES),\
		$(eval $(call extract-unmapped-pairs,$sample)))
