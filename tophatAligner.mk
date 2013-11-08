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

BAM_SUFFIX := tophat.sorted.filtered
ifeq ($(NO_REALN),false)
BAM_SUFFIX := $(BAM_SUFFIX).realn
endif

ifeq ($(DUP_TYPE),rmdup)
BAM_SUFFIX := $(BAM_SUFFIX).rmdup
else ifeq ($(DUP_TYPE),markdup) 
BAM_SUFFIX := $(BAM_SUFFIX).markdup
endif

ifeq ($(NO_RECAL),false)
BAM_SUFFIX := $(BAM_SUFFIX).recal
endif

BAM_SUFFIX := $(BAM_SUFFIX).bam

TOPHAT_BAMS = $(foreach sample,$(SAMPLES),tophat/bam/$(sample).$(BAM_SUFFIX))
tophat_bams : $(addsuffix .md5,$(TOPHAT_BAMS)) $(addsuffix .bai,$(TOPHAT_BAMS))
#tophat_unmapped_bams : $(foreach sample,$(SAMPLES),tophat/unmapped_bam/$(sample).unmapped.bam)
#tophat_junctions: $(foreach sample,$(SAMPLES),tophat/junctions/$(sample)_junctions.bed)

tophat/bam/%.tophat.bam.md5 : fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5
	$(call LSCRIPT_PARALLEL_MEM,4,6G,10G,"$(CHECK_MD5) $(TOPHAT) $(TOPHAT_OPTS) -o $(@D)/$* $(BOWTIE_REF) $(word 1,$^) $(word 2,$^) && ln -f tophat/$*/accepted_hits.bam $(@M) && $(MD5)")


include ~/share/modules/fastq.mk
include ~/share/modules/processBam.mk
