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
STAR_CHIMERIC_OPTS = --twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 \
					 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 \
					 --chimSegmentReadGapMax parameter 3 --alignSJstitchMismatchNmax 5 -1 5 5 \
					 --chimOutType WithinBAM
STAR = STAR
STAR_OPTS = --genomeDir $(STAR_REF_DIR) --outSAMtype BAM SortedByCoordinate \
			$(if $(findstring true,$(STAR_CHIMERIC)),$(STAR_CHIMERIC_OPTS)) \
			--quantMode GeneCounts


$(if $(STAR_REF_DIR),,$(error no STAR genome provided))

PHONY += star_bams star star_chimeric_junctions
star : star_bams $(if $(findstring true,$(STAR_CHIMERIC)),star_chimeric_junctions)
star_bams : $(foreach sample,$(SAMPLES),bam/$(sample).bam bam/$(sample).bam.bai)
star_chimeric_junctions : $(foreach sample,$(SAMPLES),star/$(sample).Chimeric.out.junction)

define align-split-fastq
star/$2.star_align_timestamp : $3
	$$(call RUN,-n 4 -s 6G -m 10G,"$$(STAR) $$(STAR_OPTS) \
		--outFileNamePrefix star/$2. --runThreadN 4 \
		--outSAMattrRGline \"ID:$2\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
		--readFilesIn $$^ --readFilesCommand zcat && touch $$@")
star/$2.Aligned.sortedByCoord.out.bam : star/$2.star_align_timestamp
star/$2.Chimeric.out.junction : star/$2.star_align_timestamp
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

star/%.Aligned.sortedByCoord.out.bam star/%.Chimeric.out.junction : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-n 4 -s 6G -m 10G,"$(STAR) $(STAR_OPTS) \
		--outFileNamePrefix star/$*. --runThreadN 4 \
		--outSAMattrRGline \"ID:$*\" \"LB:$*\" \"SM:$*\" \"PL:${SEQ_PLATFORM}\" \
		--readFilesIn $^ --readFilesCommand zcat")

star/bam/%.star.sorted.bam : star/%.Aligned.sortedByCoord.out.bam
	$(INIT) mv $< $@

bam/%.bam : star/bam/%.star.$(BAM_SUFFIX)
	$(INIT) ln -f $(<) $(@) 

define merged-chimeric-junction
star/$1.Chimeric.out.junction : $$(foreach split,$$(split.$1),star/$$(split).Chimeric.out.junction)
	$$(INIT) sort -V $$^ > $$@
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call merged-chimeric-junction,$(sample))))

include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
include modules/aligners/align.mk

.PHONY: $(PHONY)
.SECONDARY: 
.DELETE_ON_ERROR:
