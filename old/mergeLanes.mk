# Merge lane samples using a LANE_SAMPLE_FILE
include ~/share/modules/Makefile.inc

#SAMPLE_FILE = samples.txt
#SAMPLES = $(shell cat $(SAMPLE_FILE))

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LANE_SAMPLE_FILE ?= samples.lane.txt # generate using scripts/findLanes.sh
SAMPLES = $(shell cut -f1 $(LANE_SAMPLE_FILE) | sort | uniq)
LANES = $(shell cut -f2 $(LANE_SAMPLE_FILE))

$(foreach i,$(shell seq 1 $(words $(LANES))),$(eval lane_lookup.$(word $i,$(LANES)) := $(word $i,$(LANES))))
$(foreach i,$(shell seq 1 $(words $(LANES))),$(eval sample_lookup.$(lane_lookup.$(word $i,$(LANES))) += $(word $i,$(LANES))))

SAM_TO_FASTQ = $(JAVA) -Xmx8G -jar $(JARDIR)/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT
MERGE_BAM = $(JAVA) -Xmx8G -jar $(JARDIR)/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT

LOGDIR = log

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : merged_bam merged_fastq

merged_bam : $(foreach sample,$(SAMPLES),bam/$(sample).bam)

merged_fastq : $(foreach i,1 2,$(foreach sample,$(SAMPLES),fastq/$(sample).$i.fastq.gz))

define bam-header
merged_bam/$1.header.sam : $$(foreach lane,$2,bam/$$(lane).bam)
	$$(INIT) $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@; \
	for bam in $$^; do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@; done
endef
$(foreach sample,$(SAMPLES),$(eval $(call bam-header,$(sample),$(sample_lookup.$(sample)))))

define merged-bam
merged_bam/$1.bam : merged_bam/$1.header.sam $$(foreach lane,$2,bam/$$(lane).bam)
	$$(call INIT_MEM,12G,15G) $$(SAMTOOLS) merge -h $$< $$@ $$(filter %.bam,$$^) &> $$(LOGDIR)/$$(@F).merge.log
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(sample_lookup.$(sample)))))

define merged-fastq
fastq/$1.%.fastq.gz : $$(foreach lane,$2,fastq/$$(lane).%.fastq.gz)
	$$(INIT) zcat $^ | gzip > $@ 2> $$(LOG)
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-fastq,$(sample),$(sample_lookup.$(sample)))))


bam/%.bam : merge_bam/%.rg.rmdup.bam
	$(MKDIR) $(@D); ln -f $< $@

include ~/share/modules/processBam.mk
include ~/share/modules/fastq.mk
