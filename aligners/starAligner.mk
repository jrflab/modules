# This module is for the Tophat aligner
# input: $(SAMPLES) 
# Options: 
# 	BAM_PHRED64 = true/false
# 	NO_NOVEL_SPLICING = true/false
# 	NUM_CORES = 4
#	INNER_MATE_DIST = 200
include modules/Makefile.inc

BAM_NO_REALN = true
BAM_NO_RECAL = true
BAM_NO_FILTER = true
BAM_DUP_TYPE = none

ALIGNER := star
include modules/aligners/align.inc

LOGDIR = log/star.$(NOW)

SEQ_PLATFROM = illumina

STAR_CHIMERIC ?= true
STAR = $(HOME)/share/usr/bin/STAR
STAR_OPTS = --genomeDir $(STAR_REF_DIR) --outSAMtype BAM SortedByCoordinate \
			$(if $(findstring true,$(STAR_CHIMERIC)),--chimSegmentMin 20 --chimOutType WithinBAM) \
			--quantMode GeneCounts


$(if $(STAR_REF_DIR),,$(error no STAR genome provided))

PHONY += star_bams
star_bams : $(foreach sample,$(SAMPLES),bam/$(sample).bam bam/$(sample).bam.bai)

define align-split-fastq
star/$2.Aligned.sortedByCoord.out.bam : $3
	$$(call LSCRIPT_PARALLEL_MEM,4,6G,10G,"$$(STAR) $$(STAR_OPTS) \
		--outFileNamePrefix star/$2. --runThreadN 4 \
		--outSAMattrRGline \"ID:$2\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
		--readFilesIn $$^ --readFilesCommand zcat")
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))


star/%.Aligned.sortedByCoord.out.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,4,6G,10G,"$(STAR) $(STAR_OPTS) \
		--outFileNamePrefix star/$*. --runThreadN 4 \
		--readFilesIn $^ --readFilesCommand zcat")

star/bam/%.star.sorted.bam : star/%.Aligned.sortedByCoord.out.bam
	$(INIT) mv $< $@

bam/%.bam : star/bam/%.star.$(BAM_SUFFIX)
	$(INIT) ln -f $(<) $(@) 


include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
include modules/aligners/align.mk

.PHONY: $(PHONY)
.SECONDARY: 
.DELETE_ON_ERROR:
