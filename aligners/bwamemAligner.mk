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

..DUMMY := $(shell mkdir -p version; $(BWA) &> version/bwamem.txt; echo "options: $(BWA_ALN_OPTS)" >> version/bwamem.txt )
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bwamem


BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

bwamem : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam : bwamem/bam/%.bwamem.$(BAM_SUFFIX)
	$(call LSCRIPT,"ln -f $(<) $(@)")

ifdef SPLIT_SAMPLES
define bam-header
bwamem/sam/$1.header.sam : $$(foreach split,$2,bwamem/bam/$$(split).bwamem.sorted.bam)
	$$(INIT) $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$^; do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp
endef
$(foreach sample,$(SAMPLES),$(eval $(call bam-header,$(sample),$(split.$(sample)))))

define merged-bam
bwamem/bam/$1.bwamem.sorted.bam : bwamem/sam/$1.header.sam $$(foreach split,$2,bwamem/bam/$$(split).bwamem.sorted.bam)
	if [ `echo "$$(filter %.bam,$$(^))" | wc -w` -gt 1 ]; then \
		$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@) $$(filter %.bam,$$^) && $$(RM) $$^"); \
	else \
		ln -f $$(word 2,$$(^)) $$(@); \
	fi
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))

#$(call align-split-fastq,name,split-name,fastqs)
define align-split-fastq
bwamem/bam/$2.bwamem.bam : $3
	$$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"$$(BWA) mem -t 8 $$(BWA_ALN_OPTS) -R \"@RG\tID:$2\tLB:$1\tPL:$${SEQ_PLATFORM}\tSM:$1\" $$(REF_FASTA) $$^ | $$(SAMTOOLS) view -bhS - > $$@")
endef
$(foreach ss,$(SPLIT_SAMPLES),$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss)))))

else

bwamem/bam/%.bwamem.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"$(BWA) mem -t 8 $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

bwamem/bam/%.bwamem.bam : fastq/%.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"$(BWA) mem -t 8 $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

endif

fastq/%.fastq.gz : fastq/%.fastq
	$(call LSCRIPT,"gzip -c $< > $(@) && $(RM) $<")

include modules/bam_tools/processBam.mk
include modules/fastq_tools/fastq.mk
