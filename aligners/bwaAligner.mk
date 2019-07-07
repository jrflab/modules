include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/aligners/align.inc

ALIGNER := bwa
LOGDIR ?= log/bwa.$(NOW)

SAMTOOLS_SORT_MEM = 2000000000
SEQ_PLATFORM = illumina

VPATH ?= unprocessed_bam

# use fastq; otherwise use bams
FASTQ_CHUNKS := 10
FASTQ_CHUNK_SEQ := $(shell seq 1 $(FASTQ_CHUNKS))
FASTQUTILS = $(HOME)/share/usr/ngsutils/bin/fastqutils

BWA_ALN_OPTS ?= 
#BWA_ALN_OPTS ?= -q 20

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bwa
..DUMMY := $(shell mkdir -p version; $(BWA) &> version/bwa.txt; echo "options: $(BWA_ALN_OPTS)" >> version/bwa.txt )


BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
bwa : $(addsuffix ,$(BWA_BAMS)) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam : bwa/bam/%.bwa.$(BAM_SUFFIX)
	$(call RUN,,"ln -f $(<) $(@) ")

ifdef SPLIT_SAMPLES

.SECONDEXPANSION:
define sai-split-fastq-pair
bwa/sai/$1.$3.sai : $2
	$$(call RUN,-n 8 -s 1G -m 1.2G,"$$(BWA) aln $$(BWA_ALN_OPTS) -t 8 $$(REF_FASTA) $$(<) > $$(@)")
endef
$(foreach ss,$(SPLIT_SAMPLES), \
	$(if $(fq.$(ss)),\
	$(foreach i,1 2,\
	$(eval $(call sai-split-fastq-pair,$(ss),$(word $i,$(fq.$(ss))),$i)))))

define align-split-fastq
bwa/bam/$2.bwa.bam : bwa/sai/$2.1.sai bwa/sai/$2.2.sai $3
	$$(call RUN,-s 4G -m 10G,"$$(BWA) sampe -P -r \"@RG\tID:$2\tLB:$1\tPL:$${SEQ_PLATFORM}\tSM:$1\" $$(REF_FASTA) $$^ | $$(SAMTOOLS) view -bhS - > $$@")
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
endif

bwa/sai/%.sai : fastq/%.fastq.gz
	$(call RUN,-n 8 -s 1G -m 1.2G,"$(BWA) aln $(BWA_ALN_OPTS) -t 8 $(REF_FASTA) $(<) > $(@) ")

bwa/bam/%.bwa.bam : bwa/sai/%.1.sai bwa/sai/%.2.sai fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call RUN,-s 4G -m 10G,"$(BWA) sampe -P -r \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $(^) | $(SAMTOOLS) view -bhS - > $(@) ")

fastq/%.fastq.gz : fastq/%.fastq
	$(call RUN,,"gzip -c $< > $(@) && $(RM) $< ")


include modules/bam_tools/processBam.mk
include modules/fastq_tools/fastq.mk
include modules/aligners/align.mk
