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
EXTRACT_FASTQ ?= true
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
all : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam : bwa/bam/%.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@

fastq/%.fastq.gz : fastq/%.fastq
	$(INIT) gzip -c $< > $@ 2> $(LOG) && $(RM) $< 

#ifeq ($(SPLIT_FASTQ),true)
#$(foreach i,$(FASTQ_CHUNK_SEQ),split_fastq/%.$i.fastq.gz) : fastq/%.fastq.gz
#$(INIT) $(FASTQUTILS) split -gz $< split_fastq/$* $(FASTQ_CHUNKS) &> $(LOG)
#bwa/sai/%.sai : split_fastq/%.fastq.gz
#$(call LAUNCH_MEM,8 1G 1.2G) $(BWA) aln -t 8 $(REF_FASTA) $< > $@ 2> $(LOGDIR)/$@.log
#define bwa-chunk
#bwa/split_bam/%.$1.bwa.bam : bwa/sai/%.1.$1.sai bwa/sai/%.2.$1.sai split_fastq/%.1.$1.fastq.gz split_fastq/%.2.$1.fastq.gz
#$$(call INIT_MEM,4G,10G) \
#LBID=`echo "$$*" | sed 's/_[0-9]\+//'`; \
#$$(BWA) sampe -P -r "@RG\tID:$$*\tLB:$$$${LBID}\tPL:$${SEQ_PLATFORM}\tSM:$$$${LBID}" $$(REF_FASTA) $$^ 2> $$(LOGDIR)/$$@.log | $$(SAMTOOLS) view -bhS - > $$@ && $$(RM) $$^
#endef
#$(foreach i,$(FASTQ_CHUNK_SEQ),$(eval $(call bwa-chunk,$i)))
#bwa/bam/%.bwa.bam : $(foreach i,$(FASTQ_CHUNK_SEQ),bwa/split_bam/%.$i.bwa.sorted.bam)
#$(call INIT_MEM,2G,4G) $(SAMTOOLS) merge $@ $^ && $(RM) $^ &> $(LOG)
#else

ifeq ($(EXTRACT_FASTQ),true)
bwa/sai/%.sai : fastq/%.fastq.gz
	$(INIT) $(call LAUNCH_PARALLEL_MEM,8,1G,1.2G) "$(BWA) aln -t 8 $(REF_FASTA) $< > $@ 2> $(LOGDIR)/$@.log"
bwa/bam/%.bwa.bam : bwa/sai/%.1.sai bwa/sai/%.2.sai fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(INIT) LBID=`echo "$*" | sed 's/_[0-9]\+//'`; \
	$(call LAUNCH_MEM,4G,10G) "$(BWA) sampe -P -r \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $^ 2> $(LOGDIR)/$@.log | $(SAMTOOLS) view -bhS - > $@"
else
bwa/sai/%.1.sai : unprocessed_bam/%.bam
	$(INIT) $(call LAUNCH_PARALLEL_MEM,8,1G,1.2G) $(BWA) aln -t 8 $(REF_FASTA) -b -1 $< > $@ 2> $(LOGDIR)/$@.log
bwa/sai/%.2.sai : unprocessed_bam/%.bam
	$(INIT) $(call LAUNCH_PARALLEL_MEM,8,1G,1.2G) $(BWA) aln -t 8 $(REF_FASTA) -b -2 $< > $@ 2> $(LOGDIR)/$@.log
bwa/bam/%.bwa.bam : bwa/sai/%.1.sai bwa/sai/%.2.sai unprocessed_bam/%.bam
	$(INIT) LBID=`echo "$*" | sed 's/_[0-9]\+//'`; \
	$(call LAUNCH_MEM,4G,8G) "$(BWA) sampe -P -r \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $(word 1,$^) $(word 2, $^) $(word 3 $^) $(word 3 $^) 2> $(LOGDIR)/$@.log | $(SAMTOOLS) view -bhS - > $@" && $(RM) $<
endif



include ~/share/modules/processBam.mk
include ~/share/modules/fastq.mk
