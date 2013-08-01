# vim: set ft=make :
# BWA alignment of short reads
# OPTIONS: NO_MARKDUP = true/false (default: false)
# 		   EXTRACT_FASTQ = true/false (default: false)
# 		   NO_RECAL = true/false (default: false)

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

SAMPLE_FILE ?= samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR := log/bwa.$(NOW)

SAMTOOLS_SORT_MEM = 2000000000
SEQ_PLATFORM = illumina

VPATH ?= unprocessed_bam

# use fastq; otherwise use bams
DUP_TYPE ?= rmdup
NO_RECAL ?= false
NO_REALN ?= false
SPLIT_CHR ?= true
SPLIT_FASTQ ?= false

FASTQ_CHUNKS := 10
FASTQ_CHUNK_SEQ := $(shell seq 1 $(FASTQ_CHUNKS))
FASTQUTILS = $(HOME)/share/usr/ngsutils/bin/fastqutils

.SECONDARY:
.DELETE_ON_ERROR: 

ifeq ($(SPLIT_FASTQ),true)
BAM_SUFFIX := bwa.filtered
else
BAM_SUFFIX := bwa.sorted.filtered
endif

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

BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
all : $(addsuffix .md5,$(BWA_BAMS)) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam.md5 : bwa/bam/%.$(BAM_SUFFIX).md5
	$(INIT) cp $< $@ && ln -f $(<:.md5=) $(@:.md5=)

fastq/%.fastq.gz.md5 : fastq/%.fastq
	$(INIT) gzip -c $< > $(@:.md5=) 2> $(LOG) && $(RM) $< && $(MD5)

bwa/sai/%.sai.md5 : fastq/%.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,8,1G,1.2G,"$(BWA) aln -t 8 $(REF_FASTA) $(<:.md5=) > $(@:.md5=) 2> $(LOG) && $(MD5)")
#echo "$(BWA) aln -t 8 $(REF_FASTA) $(<:.md5=) > $(@:.md5=) 2> $(LOG) && $(MD5)" | $(call LSCRIPT_PARALLEL_MEM,8,1G,1.2G)

bwa/bam/%.bwa.bam.md5 : bwa/sai/%.1.sai.md5 bwa/sai/%.2.sai.md5 fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[0-9]\+//'`; \
	$(call LSCRIPT_MEM,4G,10G,"$(CHECK_MD5) $(BWA) sampe -P -r \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $(^:.md5=) 2> $(LOG) | $(SAMTOOLS) view -bhS - > $(@:.md5=) && $(MD5)")

include ~/share/modules/processBamMD5.mk
include ~/share/modules/fastq.mk
