# add in new fastq to existing bams
include ~/share/modules/Makefile.inc

LOGDIR ?= log/merge_fastq.$(NOW)

.PHONY: merge_fastq
.DELETE_ON_ERROR:
.SECONDARY:
.SECONDEXPANSION: 

ALIGNER ?= bwamem

merge_fastq : $(foreach sample,$(SAMPLES),merged_bam/$(sample).bam.md5)

include ~/share/modules/aligners/$(ALIGNER)Aligner.mk

# link existing bam
merged_bam/%.2.bam.md5 : bam/%.bam.md5 
	$(INIT) cp $< $@ && ln -f $(<:.md5=) $(@:.md5=)

merged_bam/%.1.bam.md5 : $(ALIGNER)/%.$(ALIGNER).$(BAM_SUFFIX).md5
	$(INIT) cp $< $@ && ln -f $(<:.md5=) $(@:.md5=)

merged_bam/%.header.sam : merged_bam/%.1.bam.md5 merged_bam/%.2.bam.md5
	$(INIT) $(SAMTOOLS) view -H $(<M) | grep -v '^@RG' > $@.tmp; \
	for bam in $(^M); do $(SAMTOOLS) view -H $$bam | grep '^@RG' >> $@.tmp; done; \
	uniq $@.tmp > $@ && $(RM) $@.tmp

merged_bam/%.bam.md5 : merged_bam/%.header.sam merged_bam/%.1.bam.md5 merged_bam/%.2.bam.md5
	$(call LSCRIPT_MEM,12G,15G,"$(SAMTOOLS) merge -f -h $< $(@M) $(filter %.bam,$(^M)) && $(MD5)")


