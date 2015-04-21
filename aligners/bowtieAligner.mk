# performs bowtie alignment from reads extracted from bam files
# INPUT: bam files
# OUTPUT: bowtie aligned bam files
# OPTIONS: PHRED64 = true/false
# 		   LOCAL = true/false (preform local alignments)
# 		   RMDUP = true/false
include modules/Makefile.inc
include modules/aligners/align.inc

VPATH ?= unprocessed_bam

LOGDIR = log/bowtie.$(NOW)

BOWTIE_OPTS = -x $(BOWTIE_REF) 

LOCAL ?= false
PHRED64 ?= false
SEQ_PLATFORM ?= ILLUMINA
NUM_CORES ?= 4

ifeq ($(PHRED64),true)
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

bowtie_bams : $(addsuffix .md5,$(BAMS)) $(addsuffix .bai,$(BAMS))

# memory for human genome: ~3.2G
bowtie/bam/%.bwt.bam.md5 : fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5
	$(call LSCRIPT_PARALLEL_MEM,4,1G,1.5G,"$(CHECK_MD5) LBID=`echo \"$*\" | sed 's/_[A-Za-z0-9]\+//'`; \
		$(BOWTIE) $(BOWTIE_OPTS) --rg-id $* --rg \"LB:\$${LBID}\" --rg \"PL:${SEQ_PLATFORM}\" --rg \"SM:\$${LBID}\" -p $(NUM_CORES) --1 $(word 1,$(^:.md5=)) -2 $(word 2,$(^:.md5=)) | $(SAMTOOLS) view -bhS - > $(@:.md5=) && $(MD5)")

bowtie/bam/%.bwt.bam.md5 : fastq/%.fastq.gz.md5
	$(call LSCRIPT_PARALLEL_MEM,4,1G,1.5G,"$(CHECK_MD5) LBID=`echo \"$*\" | sed 's/_[A-Za-z0-9]\+//'`; \
		$(BOWTIE) $(BOWTIE_OPTS) --rg-id $* --rg \"LB:\$${LBID}\" --rg \"PL:${SEQ_PLATFORM}\" --rg \"SM:\$${LBID}\" -p $(NUM_CORES) -U $(<:.md5=) | $(SAMTOOLS) view -bhS - > $(@:.md5=) && $(MD5)")


bam/%.bam.md5 : bowtie/bam/%.bwt.$(BAM_SUFFIX).md5
	$(call LSCRIPT,"ln -f $(<:.md5=) $(@:.md5=) && $(MD5)")

ifdef SPLIT_SAMPLES
define merged-bam
ifeq ($(shell echo "$(words $2) > 1" | bc),1)
bowtie/bam/$1.header.sam : $$(foreach split,$2,bowtie/bam/$$(split).bwt.sorted.bam.md5)
	$$(INIT) $$(SAMTOOLS) view -H $$(<M) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp

bowtie/bam/$1.bwt.sorted.bam.md5 : bowtie/bam/$1.header.sam $$(foreach split,$2,bowtie/bam/$$(split).bwt.sorted.bam.md5)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@M) $$(filter %.bam,$$(^M)) && $$(MD5) && $$(RM) $$(^M) $$^")
endif
ifeq ($(shell echo "$(words $2) == 1" | bc),1)
bowtie/bam/$1.bwt.bam.md5 : bowtie/bam/$2.bwt.bam.md5
	$$(INIT) mv $$(<M) $$(@M) && $$(MD5)
endif
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))
endif


include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
