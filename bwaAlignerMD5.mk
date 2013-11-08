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

BWA_ALN_OPTS ?= 
#BWA_ALN_OPTS ?= -q 20

.SECONDARY:
.DELETE_ON_ERROR: 

BAM_SUFFIX := bwa.sorted.filtered

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

ifdef SPLIT_SAMPLES
define bam-header
bwa/bam/$1.header.sam : $$(foreach split,$2,bwa/bam/$$(split).bwa.sorted.bam.md5)
	$$(INIT) $$(SAMTOOLS) view -H $$(<M) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split_lookup.$(sample)))))

define merged-bam
bwa/bam/$1.bwa.sorted.bam.md5 : bwa/bam/$1.header.sam $$(foreach split,$2,bwa/bam/$$(split).bwa.sorted.bam.md5)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@M) $$(filter %.bam,$$(^M)) && $$(MD5) && $$(RM) $$(^M) $$^")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))
endif

fastq/%.fastq.gz.md5 : fastq/%.fastq
	$(call LSCRIPT,"gzip -c $< > $(@:.md5=) && $(RM) $< && $(MD5)")

bwa/sai/%.sai.md5 : fastq/%.fastq.gz.md5
	$(call LSCRIPT_PARALLEL_MEM,8,1G,1.2G,"$(CHECK_MD5) $(BWA) aln $(BWA_ALN_OPTS) -t 8 $(REF_FASTA) $(<:.md5=) > $(@:.md5=) 2> $(LOG) && $(MD5)")
#echo "$(BWA) aln -t 8 $(REF_FASTA) $(<:.md5=) > $(@:.md5=) 2> $(LOG) && $(MD5)" | $(call LSCRIPT_PARALLEL_MEM,8,1G,1.2G)

bwa/bam/%.bwa.bam.md5 : bwa/sai/%.1.sai.md5 bwa/sai/%.2.sai.md5 fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5
	LBID=`echo "$*" | sed 's/_[0-9]\+//'`; \
	$(call LSCRIPT_MEM,4G,10G,"$(CHECK_MD5) $(BWA) sampe -P -r \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $(^:.md5=) 2> $(LOG) | $(SAMTOOLS) view -bhS - > $(@:.md5=) && $(MD5)")

include ~/share/modules/processBamMD5.mk
include ~/share/modules/fastq.mk
