# Merge samples using SAMPLE_SPLIT_FILE entries and then fix the read groups
include modules/Makefile.inc

LOGDIR = log/merge.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : merged_bam

merged_bam : $(foreach sample,$(SPLIT_SAMPLES),bam/$(sample).bam bam/$(sample).bam.bai)

define merged-bam
merged_bam/$1.header.sam : $$(foreach split,$$(split_lookup.$1),bam/$$(split).bam)
	$$(INIT) $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$^; do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp

merged_bam/$1.bam : merged_bam/$1.header.sam $$(foreach split,$$(split_lookup.$1),bam/$$(split).bam)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@) $$(filter %.bam,$$^)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample))))

bam/%.bam : merged_bam/%.rg.bam
	$(INIT) ln -f $< $@

include modules/bam_tools/processBam.mk
