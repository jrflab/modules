# add in new fastq to existing bams
include modules/Makefile.inc

LOGDIR ?= log/merge_fastq.$(NOW)

.PHONY: merge_fastq
.DELETE_ON_ERROR:
.SECONDARY:

MERGE_SAMPLE_FILE ?= merge_samples.txt
ifneq ($(wildcard $(MERGE_SAMPLE_FILE)),)
  MERGE_SAMPLES ?= $(shell sed '/^\#/d' $(MERGE_SAMPLE_FILE))
endif

ALIGNER ?= bwamem

merge_fastq : $(foreach sample,$(MERGE_SAMPLES),$(if "$(wildcard bam/$(sample).bam)",merged_bam/$(sample).bam.md5,bam/$(sample).bam.md5))
	$(INIT) for bamMd5 in $(filter merged_bam/%.bam.md5,$^); do \
		cp $$bamMd5 bam/$$(basename $$bamMd5) && \
		ln -f $${bamMd5%%.md5} bam/$$(basename $${bamMd5%%.md5}); \
		done

include modules/aligners/$(ALIGNER)Aligner.mk

merged_bam/%.1.bam.md5 merged_bam/%.2.bam.md5 : $(ALIGNER)/bam/%.$(ALIGNER).$(BAM_SUFFIX).md5
	$(INIT) cp $< merged_bam/$*.1.bam.md5 && ln -f $(<M) merged_bam/$*.1.bam && \
	cp bam/$*.bam.md5 merged_bam/$*.2.bam.md5 && ln -f bam/$*.bam merged_bam/$*.2.bam

merged_bam/%.header.sam : merged_bam/%.1.bam.md5 merged_bam/%.2.bam.md5
	$(INIT) $(SAMTOOLS) view -H $(<M) | grep -v '^@RG' > $@.tmp; \
	for bam in $(^M); do $(SAMTOOLS) view -H $$bam | grep '^@RG' >> $@.tmp; done; \
	uniq $@.tmp > $@ && $(RM) $@.tmp

merged_bam/%.bam.md5 : merged_bam/%.header.sam merged_bam/%.1.bam.md5 merged_bam/%.2.bam.md5
	$(call LSCRIPT_MEM,12G,15G,"$(SAMTOOLS) merge -f -h $< $(@M) $(filter %.bam,$(^M)) && $(MD5)")


