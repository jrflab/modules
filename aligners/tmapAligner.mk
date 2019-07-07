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

bam/%.bam : tmap/bam/%.tmap.$(BAM_SUFFIX)
	$(INIT) cp $< $@ && ln -f $(<) $(@)

tmap/sam/%.header.sam : unprocessed_bam/%.bam
	$(INIT) $(SAMTOOLS) view -H $< | grep -e '^@HD' -e '^@RG' > $@

tmap/bam/%.$(TMAP_MODE).bam : tmap/sam/%.header.sam unprocessed_bam/%.bam
	$(call RUN,-c -n 4 -s 6G -m 8G,"$(SAMTOOLS) reheader $^ | $(TMAP) $(TMAP_MODE) $(TMAP_OPTS) -Q 2 -f $(REF_FASTA) -i bam -s $(@) -o 1 -n 4 ")

define align-split-fastq
tmap/bam/$2.tmap.bam : $3
	$$(call RUN,-c -n 4 -s 6G -m 8G,"zcat $$< | $$(TMAP) $$(TMAP_MODE) $$(TMAP_OPTS) -f $$(REF_FASTA) -i fastq -s $$(@) -Q 0 -o 1 -n 4 -R ID:$2 -R SM:$1 -R PL:$$(SEQ_PLATFORM) -R PU:00000000")
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)), \
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

tmap/bam/%.tmap.bam : fastq/%.fastq.gz
	$(call RUN,-c -n 4 -s 6G -m 8G,"zcat $< | $(TMAP) $(TMAP_MODE) $(TMAP_OPTS) -f $(REF_FASTA) -i fastq -s $(@) -Q 0 -o 1 -n 4 -R ID:$* -R SM:$* -R PL:$(SEQ_PLATFORM) -R PU:00000000 ")

include modules/bam_tools/processBam.mk
include modules/fastq_tools/fastq.mk
include modules/aligners/align.mk
