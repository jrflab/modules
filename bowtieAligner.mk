# performs bowtie alignment from reads extracted from bam files
# INPUT: bam files
# OUTPUT: bowtie aligned bam files
# OPTIONS: PHRED64 = true/false
# 		   LOCAL = true/false (preform local alignments)
# 		   RMDUP = true/false
include ~/share/modules/Makefile.inc

VPATH ?= gsc_bam
SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

LOGDIR = log/bowtie.$(NOW)

OPTS = -x $(BOWTIE_UCSC_REF) 

LOCAL ?= false
PHRED64 ?= false
RMDUP ?= true
SEQ_PLATFORM ?= ILLUMINA
NUM_CORES ?= 10

ifeq ($(PHRED64),true)
  OPTS += --phred64
endif

ifeq ($(LOCAL),true)
  OPTS += --local
endif

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

.PHONY : all bowtie_bams

all : bowtie_bams
bowtie_bams : $(BAMS) $(addsuffix .bai,$(BAMS))

# memory for human genome: ~3.2G
bowtie/bam/%.bowtie.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[0-9]\+//'`; \
	$(call INIT_PARALLEL_MEM,10,0.4G,0.6G) $(BOWTIE) $(OPTS) --rg-id $* --rg "LB:$${LBID}" --rg "PL:${SEQ_PLATFORM}" --rg "SM:$${LBID}" -p $(NUM_CORES) --1 $(word 1,$^) -2 $(word 2,$^) 2> $(LOGDIR)/$(@F).log | $(SAMTOOLS) view -bhS - > $@ 

ifeq ($(RMDUP),true) 
bam/%.bam : bowtie/bam/%.bowtie.sorted.filtered.rmdup.bam
	$(MKDIRS); ln -f $< $@
else
bam/%.bam : bowtie/bam/%.bowtie.sorted.filtered.bam
	$(MKDIRS); ln -f $< $@
endif

include ~/share/modules/fastq.mk
include ~/share/modules/processBam.mk
