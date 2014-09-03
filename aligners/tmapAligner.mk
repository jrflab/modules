# vim: set ft=make :
# TMAP alignment of short reads
# OPTIONS: NO_MARKDUP = true/false (default: false)
# 		   EXTRACT_FASTQ = true/false (default: false)
# 		   NO_RECAL = true/false (default: false)

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

LOGDIR := log/tmap.$(NOW)

SAMTOOLS_SORT_MEM = 2000000000

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

TMAP = $(HOME)/share/usr/bin/tmap
TMAP_MODE ?= map3
TMAP_OPTS = -Q 2

.SECONDARY:
.DELETE_ON_ERROR: 

BAM_SUFFIX := $(TMAP_MODE).sorted.filtered

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

TMAP_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
all : $(addsuffix .md5,$(TMAP_BAMS)) $(addsuffix .bai,$(TMAP_BAMS))

bam/%.bam.md5 : tmap/bam/%.$(BAM_SUFFIX).md5
	$(INIT) cp $< $@ && ln -f $(<:.md5=) $(@:.md5=)

tmap/sam/%.header.sam : unprocessed_bam/%.bam
	$(INIT) $(SAMTOOLS) view -H $< | grep -e '^@HD' -e '^@RG' > $@

tmap/bam/%.$(TMAP_MODE).bam.md5 : tmap/sam/%.header.sam unprocessed_bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,6G,8G,"$(SAMTOOLS) reheader $^ | $(TMAP) $(TMAP_MODE) $(TMAP_OPTS) -f $(REF_FASTA) -i bam -s $(@M) -o 1 -n 4 && $(MD5)")

ifdef SPLIT_SAMPLES
define bam-header
tmap/bam/$1.header.sam : $$(foreach split,$2,tmap/bam/$$(split).$$(TMAP_MODE).sorted.bam.md5)
	$$(INIT) $$(SAMTOOLS) view -H $$(<M) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split_lookup.$(sample)))))

define merged-bam
tmap/bam/$1.$$(TMAP_MODE).sorted.bam.md5 : tmap/bam/$1.header.sam $$(foreach split,$2,tmap/bam/$$(split).$$(TMAP_MODE).sorted.bam.md5)
	if [ `echo "$$(filter %.bam,$$(^M))" | wc -w` -gt 1 ]; then \
		$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@M) $$(filter %.bam,$$(^M)) && $$(MD5) && $$(RM) $$(^M) $$^"); \
	else \
		ln -f $$(word 2,$$(^M)) $$(@M) && ln -f $$(word 2,$$^) $$@; \
	fi
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))
endif

include ~/share/modules/processBam.mk
include ~/share/modules/fastq.mk
