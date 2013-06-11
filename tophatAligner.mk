# This module is for the Tophat aligner
# input: $(SAMPLES) 
# Options: 
# 	PHRED64 = true/false
# 	NO_NOVEL_SPLICING = true/false
# 	NUM_CORES = 4
#	INNER_MATE_DIST = 200
include ~/share/modules/Makefile.inc

SAMPLE_FILE ?= samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

NUM_CORES ?= 4
INNER_MATE_DIST ?= 200
NO_NOVEL_SPLICING ?= false 

TOPHAT_OPTS = --mate-inner-dist $(INNER_MATE_DIST) -G $(REFSEQ_GTF) -p ${NUM_CORES}
#TOPHAT_SGE_RREQ = $(call MEM_FREE,4G,5G) -q all.q -pe $(PARALLEL_ENV) $(NUM_CORES) -now n
TOPHAT_SGE_RREQ = $(call MEM_FREE,6G,8G) -q all.q -pe $(PARALLEL_ENV) $(NUM_CORES) -now n

ifeq ($(PHRED64),true)
	TOPHAT_OPTS += --solexa1.3-quals
endif
ifeq ($(NO_NOVEL_SPLICING),true)
	TOPHAT_OPTS += --no-novel-juncs
endif

.PHONY : all tophat_bams tophat_junctions
	
all : tophat_bams tophat_junctions
tophat_bams : $(foreach sample,$(SAMPLES),tophat/bam/$(sample).bam) $(foreach sample,$(SAMPLES),tophat/bam/$(sample).bam.bai)
tophat_unmapped_bams : $(foreach sample,$(SAMPLES),tophat/unmapped_bam/$(sample).unmapped.bam)
tophat_junctions: $(foreach sample,$(SAMPLES),tophat/junctions/$(sample)_junctions.bed)

tophat/%.timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	SGE_RREQ="-N X$(@F) $(TOPHAT_SGE_RREQ)" \
	$(MKDIR) $(@D)/logs;\
	$(TOPHAT) $(TOPHAT_OPTS) -o $(@D)/$* $(BOWTIE_REF) $(word 1,$^) $(word 2,$^) &> $(@D)/logs/$(@F).log && touch $@

tophat/%/accepted_hits.bam : tophat/%.timestamp

tophat/unmapped_bam/%.unmapped.bam : tophat/%/unmapped.bam
	$(MKDIR) $(@D); ln -f $< $@

tophat/bam/%.bam : tophat/%/accepted_hits.sorted.filtered.rmdup.bam
	$(MKDIR) $(@D); ln -f $< $@

tophat/junctions/%_junctions.bed : tophat/%.timestamp
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,100M,1G)" \
	$(MKDIR) $(@D)/logs;\
	ln -f tophat/$*/junctions.bed $@; touch $@

include ~/share/modules/fastq.mk
include ~/share/modules/processBam.mk
