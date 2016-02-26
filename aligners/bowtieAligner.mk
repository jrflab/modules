# performs bowtie alignment from reads extracted from bam files
# INPUT: bam files
# OUTPUT: bowtie aligned bam files
# OPTIONS: BAM_PHRED64 = true/false
# 		   LOCAL = true/false (preform local alignments)
# 		   RMDUP = true/false
include modules/Makefile.inc
include modules/aligners/align.inc

VPATH ?= unprocessed_bam

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

ifdef SPLIT_SAMPLES
define merged-bam
ifeq ($(shell echo "$(words $2) > 1" | bc),1)
bowtie/bam/$1.header.sam : $$(foreach split,$2,bowtie/bam/$$(split).bwt.sorted.bam)
	$$(INIT) $$(SAMTOOLS) view -H $$(<M) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp

bowtie/bam/$1.bwt.sorted.bam : bowtie/bam/$1.header.sam $$(foreach split,$2,bowtie/bam/$$(split).bwt.sorted.bam)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@) $$(filter %.bam,$$(^)) && $$(RM) $$(^)")
endif
ifeq ($(shell echo "$(words $2) == 1" | bc),1)
bowtie/bam/$1.bwt.bam : bowtie/bam/$2.bwt.bam
	$$(INIT) mv $$(<) $$(@) 
endif
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))

define align-split-fastq
bowtie/bam/$2.bwt.bam : $3
	$$(call LSCRIPT_PARALLEL_MEM,4,1G,1.5G,"$$(BOWTIE) $$(BOWTIE_OPTS) \
		--rg-id $2 --rg \"LB:$1\" --rg \"PL:$${SEQ_PLATFORM}\" --rg \"SM:$1\" \
		-p $$(NUM_CORES) $$(if $$(<<),-1 $$(<) -2 $$(<<),-U $$<) | $$(SAMTOOLS) view -bhS - > $$(@)")
endef
$(foreach ss,$(SPLIT_SAMPLES),$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss)))))
endif


include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
