# vim: set ft=make :
# TMAP alignment of short reads
# OPTIONS: NO_MARKDUP = true/false (default: false)
# 		   EXTRACT_FASTQ = true/false (default: false)
# 		   BAM_NO_RECAL = true/false (default: false)

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/aligners/align.inc

ALIGNER := tmap
LOGDIR := log/tmap.$(NOW)

SAMTOOLS_SORT_MEM = 2000000000

FASTQ_CHUNKS := 10
FASTQ_CHUNK_SEQ := $(shell seq 1 $(FASTQ_CHUNKS))
FASTQUTILS = $(HOME)/share/usr/ngsutils/bin/fastqutils

TMAP = $(HOME)/share/usr/bin/tmap
TMAP_MODE ?= map3
TMAP_OPTS =

SEQ_PLATFORM = IONTORRENT

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: tmap

TMAP_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
tmap : $(TMAP_BAMS) $(addsuffix .bai,$(TMAP_BAMS))

bam/%.bam : tmap/bam/%.$(TMAP_MODE).$(BAM_SUFFIX)
	$(INIT) cp $< $@ && ln -f $(<) $(@)

tmap/sam/%.header.sam : unprocessed_bam/%.bam
	$(INIT) $(SAMTOOLS) view -H $< | grep -e '^@HD' -e '^@RG' > $@

tmap/bam/%.$(TMAP_MODE).bam : tmap/sam/%.header.sam unprocessed_bam/%.bam
	$(call LSCRIPT_CHECK_PARALLEL_MEM,4,6G,8G,"$(SAMTOOLS) reheader $^ | $(TMAP) $(TMAP_MODE) $(TMAP_OPTS) -Q 2 -f $(REF_FASTA) -i bam -s $(@) -o 1 -n 4 ")

ifdef SPLIT_SAMPLES
define bam-header
tmap/bam/$1.header.sam : $$(foreach split,$2,tmap/bam/$$(split).$$(TMAP_MODE).sorted.bam)
	$$(INIT) $$(SAMTOOLS) view -H $$(<) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp
endef
$(foreach sample,$(SAMPLES),$(eval $(call bam-header,$(sample),$(split.$(sample)))))

define merged-bam
tmap/bam/$1.$$(TMAP_MODE).sorted.bam : tmap/bam/$1.header.sam $$(foreach split,$2,tmap/bam/$$(split).$$(TMAP_MODE).sorted.bam)
	if [ `echo "$$(filter %.bam,$$(^))" | wc -w` -gt 1 ]; then \
		$$(call LSCRIPT_CHECK_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@) $$(filter %.bam,$$(^))  && $$(RM) $$^"); \
	else \
		ln -f $$(word 2,$$(^)) $$(@) \
	fi
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))

define align-split-fastq
tmap/bam/$2.$(TMAP_MODE).bam : $3
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,6G,8G,"zcat $$< | $$(TMAP) $$(TMAP_MODE) $$(TMAP_OPTS) -f $$(REF_FASTA) -i fastq -s $$(@) -Q 0 -o 1 -n 4 -R ID:$2 -R SM:$1 -R PL:$$(SEQ_PLATFORM) -R PU:00000000")
endef
$(foreach ss,$(SPLIT_SAMPLES),$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss)))))
else
tmap/bam/%.$(TMAP_MODE).bam : fastq/%.fastq.gz
	$(call LSCRIPT_CHECK_PARALLEL_MEM,4,6G,8G,"zcat $< | $(TMAP) $(TMAP_MODE) $(TMAP_OPTS) -f $(REF_FASTA) -i fastq -s $(@) -Q 0 -o 1 -n 4 -R ID:$* -R SM:$* -R PL:$(SEQ_PLATFORM) -R PU:00000000 ")

endif

include modules/bam_tools/processBam.mk
include modules/fastq_tools/fastq.mk
