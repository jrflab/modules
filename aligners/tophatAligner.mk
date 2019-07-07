include modules/Makefile.inc

BAM_NO_REALN = true
BAM_NO_RECAL = true
BAM_NO_FILTER = true
BAM_DUP_TYPE = none

ALIGNER := tophat
include modules/aligners/align.inc

LOGDIR = log/tophat.$(NOW)

TOPHAT_NUM_CORES ?= 4
NO_NOVEL_SPLICING ?= false 

TOPHAT_OPTS = --keep-fasta-order --no-sort-bam -G $(GENES_GTF) -p ${TOPHAT_NUM_CORES} --tmp-dir $(TMPDIR)/$*

SEQ_PLATFORM ?= illumina

ifeq ($(BAM_PHRED64),true)
	TOPHAT_OPTS += --solexa1.3-quals
endif
ifeq ($(NO_NOVEL_SPLICING),true)
	TOPHAT_OPTS += --no-novel-juncs
endif

..DUMMY := $(shell mkdir -p version; $(TOPHAT) --version &> version/tophat.txt; echo "options: $(TOPHAT_OPTS)" >> version/tophat.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : tophat_bams
	
BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
tophat_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
#tophat_unmapped_bams : $(foreach sample,$(SAMPLES),tophat/unmapped_bam/$(sample).unmapped.bam)
#tophat_junctions: $(foreach sample,$(SAMPLES),tophat/junctions/$(sample)_junctions.bed)

bam/%.bam : tophat/bam/%.tophat.$(BAM_SUFFIX)
	$(INIT) ln -f $(<) $(@) 

tophat/bam/%.tophat.sorted.bam : tophat/%/accepted_hits.sorted.bam tophat/%/unmapped.sorted.bam tophat/%/accepted_hits.sorted.bam.bai tophat/%/unmapped.sorted.bam.bai
	$(call RUN,-s 7G -m 7G,"$(SAMTOOLS) merge -f $(@) $(<) $(<<)")

define align-split-fastq
tophat/$2/accepted_hits.bam : $3
	$$(call RUN,-N $2_tophat -n 4 -s 6G -m 10G,"$$(TOPHAT) \
		$$(TOPHAT_OPTS) \
		--rg-id $2 --rg-library $1 \
		--rg-platform $$(SEQ_PLATFORM) --rg-sample $1 \
		-o tophat/$2 $$(BOWTIE_REF) $$(<) $$(<<)")
tophat/$2/unmapped.bam : tophat/$2/accepted_hits.bam
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

tophat/%/accepted_hits.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-N $*_tophat -n 4 -s 6G -m 10G,"LBID=`echo \"$*\" | sed 's/_.\+//'`; \
		$(TOPHAT) $(TOPHAT_OPTS) \
		--rg-id $* --rg-library \"\$${LBID}\" \
		--rg-platform \"${SEQ_PLATFORM}\" --rg-sample \"\$${LBID}\" \
		-o tophat/$* $(BOWTIE_REF) $(<) $(<<)")

include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
include modules/aligners/align.mk
