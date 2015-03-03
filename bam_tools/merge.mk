# Merge lane samples using a LANE_SAMPLE_FILE
include modules/Makefile.inc

SAM_TO_FASTQ = $(JAVA) -Xmx8G -jar $(JARDIR)/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT
MERGE_BAM = $(JAVA) -Xmx8G -jar $(JARDIR)/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT

LOGDIR = log/merge.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : merged_bam

merged_bam : $(foreach sample,$(SPLIT_SAMPLES),bam/$(sample).bam.md5) $(foreach sample,$(SPLIT_SAMPLES),bam/$(sample).bam.bai) 

define bam-header
merged_bam/$1.header.sam : $$(foreach split,$2,bam/$$(split).bam)
	$$(INIT) $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$^; do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split_lookup.$(sample)))))

define merged-bam
bam/$1.bam.md5 : merged_bam/$1.header.sam $$(foreach split,$2,bam/$$(split).bam)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@M) $$(filter %.bam,$$^) && $$(MD5)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))

include modules/bam_tools/processBam.mk
