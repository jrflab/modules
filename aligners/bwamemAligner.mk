# vim: set ft=make :
# BWA-mem alignment of short reads
# OPTIONS: NO_MARKDUP = true/false (default: false)
# 		   EXTRACT_FASTQ = true/false (default: false)
# 		   NO_RECAL = true/false (default: false)

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc
include ~/share/modules/aligners/align.inc

ALIGNER := bwamem

LOGDIR := log/bwamem.$(NOW)

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

bwamem : $(addsuffix .md5,$(BWA_BAMS)) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam.md5 : bwamem/%.bwamem.$(BAM_SUFFIX).md5
	$(INIT) cp $< $@ && ln -f $(<:.md5=) $(@:.md5=)

ifdef SPLIT_SAMPLES
define bam-header
bwamem/$1.header.sam : $$(foreach split,$2,bwamem/$$(split).bwamem.sorted.bam.md5)
	$$(INIT) $$(SAMTOOLS) view -H $$(<M) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split_lookup.$(sample)))))

define merged-bam
bwamem/$1.bwamem.sorted.bam.md5 : bwamem/$1.header.sam $$(foreach split,$2,bwamem/$$(split).bwamem.sorted.bam.md5)
	if [ `echo "$$(filter %.bam,$$(^M))" | wc -w` -gt 1 ]; then \
		$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@M) $$(filter %.bam,$$(^M)) && $$(MD5) && $$(RM) $$(^M) $$^"); \
	else \
		ln -f $$(word 2,$$(^M)) $$(@M) && ln -f $$(word 2,$$^) $$@; \
	fi
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))
endif

fastq/%.fastq.gz.md5 : fastq/%.fastq
	$(call LSCRIPT,"gzip -c $< > $(@:.md5=) && $(RM) $< && $(MD5)")

bwamem/%.bwamem.bam.md5 : fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"$(CHECK_MD5) $(BWA) mem -t 8 $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $(^M) | $(SAMTOOLS) view -bhS - > $(@:.md5=) && $(MD5)")

include ~/share/modules/bam_tools/processBam.mk
include ~/share/modules/fastq_tools/fastq.mk
