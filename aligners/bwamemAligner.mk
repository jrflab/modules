# vim: set ft=make :
# BWA-mem alignment of short reads
# OPTIONS: NO_MARKDUP = true/false (default: false)
# 		   EXTRACT_FASTQ = true/false (default: false)
# 		   BAM_NO_RECAL = true/false (default: false)

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/aligners/align.inc

ALIGNER := bwamem

LOGDIR ?= log/bwamem.$(NOW)

SAMTOOLS_SORT_MEM = 2000000000
SEQ_PLATFORM = illumina

VPATH ?= unprocessed_bam

# use fastq; otherwise use bams
FASTQ_CHUNKS := 10
FASTQ_CHUNK_SEQ := $(shell seq 1 $(FASTQ_CHUNKS))
FASTQUTILS = $(HOME)/share/usr/ngsutils/bin/fastqutils

BWA_ALN_OPTS ?= -M
#BWA_ALN_OPTS ?= -q 20
BWAMEM_REF_FASTA ?= $(REF_FASTA)

..DUMMY := $(shell mkdir -p version; $(BWA) &> version/bwamem.txt; echo "options: $(BWA_ALN_OPTS)" >> version/bwamem.txt )
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bwamem


BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

bwamem : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam : bwamem/bam/%.bwamem.$(BAM_SUFFIX)
	$(call LSCRIPT,"ln -f $(<) $(@)")

#$(call align-split-fastq,name,split-name,fastqs)
define align-split-fastq
bwamem/bam/$2.bwamem.bam : $3
	$$(call LSCRIPT_PARALLEL_MEM,8,1G,$$(if $$(findstring true,$$(PDX),4G,2G)),"$$(BWA) mem -t 8 $$(BWA_ALN_OPTS) -R \"@RG\tID:$2\tLB:$1\tPL:$${SEQ_PLATFORM}\tSM:$1\" $$(BWAMEM_REF_FASTA) $$^ | $$(SAMTOOLS) view -bhS - > $$@")
endef
$(foreach ss,$(SPLIT_SAMPLES),$(if $(fq.$(ss)),$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

bwamem/bam/%.bwamem.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"$(BWA) mem -t 8 $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(BWAMEM_REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

bwamem/bam/%.bwamem.bam : fastq/%.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"$(BWA) mem -t 8 $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(BWAMEM_REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

fastq/%.fastq.gz : fastq/%.fastq
	$(call LSCRIPT,"gzip -c $< > $(@) && $(RM) $<")

include modules/bam_tools/processBam.mk
include modules/fastq_tools/fastq.mk
include modules/aligners/align.mk
