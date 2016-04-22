# performs bowtie alignment from reads extracted from bam files
# INPUT: bam files
# OUTPUT: bowtie aligned bam files
# OPTIONS: BAM_PHRED64 = true/false
# 		   LOCAL = true/false (preform local alignments)
# 		   RMDUP = true/false
include modules/Makefile.inc
include modules/aligners/align.inc

VPATH ?= unprocessed_bam

ALIGNER := bowtie

LOGDIR = log/bowtie.$(NOW)

BOWTIE_OPTS = -x $(BOWTIE_REF) 

LOCAL ?= false
BAM_PHRED64 ?= false
SEQ_PLATFORM ?= ILLUMINA
NUM_CORES ?= 4

ifeq ($(BAM_PHRED64),true)
  BOWTIE_OPTS += --phred64
endif

ifeq ($(LOCAL),true)
  BOWTIE_OPTS += --local
endif


BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)


..DUMMY := $(shell mkdir -p version; $(BOWTIE) --version $(BOWTIE_OPTS) > version/bowtie.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: bowtie_bams

bowtie_bams : $(BAMS) $(addsuffix .bai,$(BAMS))

# memory for human genome: ~3.2G
bowtie/bam/%.bwt.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,4,1G,1.5G,"LBID=`echo \"$*\" | sed 's/_[A-Za-z0-9]\+//'`; \
		$(BOWTIE) $(BOWTIE_OPTS) --rg-id $* --rg \"LB:\$${LBID}\" --rg \"PL:${SEQ_PLATFORM}\" --rg \"SM:\$${LBID}\" -p $(NUM_CORES) \
		-1 $< -2 $(<<) | $(SAMTOOLS) view -bhS - > $(@)")

bowtie/bam/%.bwt.bam : fastq/%.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,4,1G,1.5G,"LBID=`echo \"$*\" | sed 's/_[A-Za-z0-9]\+//'`; \
		$(BOWTIE) $(BOWTIE_OPTS) --rg-id $* --rg \"LB:\$${LBID}\" --rg \"PL:${SEQ_PLATFORM}\" --rg \"SM:\$${LBID}\" -p $(NUM_CORES) \
		-U $(<) | $(SAMTOOLS) view -bhS - > $(@) ")


bam/%.bam : bowtie/bam/%.bwt.$(BAM_SUFFIX)
	$(call LSCRIPT,"ln -f $(<) $(@) ")

define align-split-fastq
bowtie/bam/$2.bwt.bam : $3
	$$(call LSCRIPT_PARALLEL_MEM,4,1G,1.5G,"$$(BOWTIE) $$(BOWTIE_OPTS) \
		--rg-id $2 --rg \"LB:$1\" --rg \"PL:$${SEQ_PLATFORM}\" --rg \"SM:$1\" \
		-p $$(NUM_CORES) $$(if $$(<<),-1 $$(<) -2 $$(<<),-U $$<) | $$(SAMTOOLS) view -bhS - > $$(@)")
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
include modules/aligners/align.mk
