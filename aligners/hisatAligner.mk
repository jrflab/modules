# This module is for the hisat aligner
# input: $(SAMPLES) 
include modules/Makefile.inc

BAM_NO_REALN = true
BAM_NO_RECAL = true
BAM_NO_FILTER = true
BAM_DUP_TYPE = none
include modules/aligners/align.inc

LOGDIR = log/hisat.$(NOW)

HISAT_NUM_CORES ?= 8
HISAT_OPTS = -q -p $(HISAT_NUM_CORES) -x $(HISAT_REF)

SEQ_PLATFORM ?= illumina
ifeq ($(SEQ_PLATFORM),true)
	HISAT_OPTS += --fr
endif

ifeq ($(PHRED64),true)
	HISAT_OPTS += --phred64
endif

..DUMMY := $(shell mkdir -p version; $(HISAT) --version &> version/hisat.txt; echo "options: $(HISAT_OPTS)" >> version/hisat.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : hisat_bams
	
BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
hisat_bams : $(BAMS) $(addsuffix .bai,$(BAMS))

bam/%.bam : hisat/bam/%.hisat.$(BAM_SUFFIX)
	$(INIT) ln -f $(<) $(@) 

hisat/bam/%.hisat.bam hisat/unmapped_bam/%.hisat_unmapped.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_NAMED_PARALLEL_MEM,$*_hisat,8,1G,1.5G,"mkdir -p hisat/bam hisat/unmapped_bam; \
		LBID=`echo \"$*\" | sed 's/_.\+//'`; \
		$(HISAT) $(HISAT_OPTS) \
		--rg-id $* --rg \"SM:\$${LBID}\" \
		--rg PL:${SEQ_PLATFORM} --rg \"LB:\$${LBID}\" \
		-S >($(SAMTOOLS2) view -bh - > hisat/bam/$*.hisat.bam) \
		--un >($(SAMTOOLS2) view -bh - > hisat/unmapped_bam/$*.hisat_unmapped.bam) \
		-1 $(<) -2 $(<<)")

ifdef SPLIT_SAMPLES
define merged-bam
hisat/bam/$1.header.sam : $$(foreach split,$2,hisat/bam/$$(split).hisat.sorted.bam)
	$$(INIT) $$(SAMTOOLS) view -H $$(<) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp

hisat/bam/$1.hisat.sorted.bam : hisat/bam/$1.header.sam $$(foreach split,$2,hisat/bam/$$(split).hisat.sorted.bam)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@) $$(filter %.bam,$$(^))  && $$(RM) $$(^)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))

# usage $(eval $(call align-merged-split-fastq,merged-sample,split,fastq(s)))
define align-merged-split-fastq
hisat/bam/$2.hisat.bam : $3
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$2_hisat,8,1G,1.5G,"mkdir -p hisat/bam hisat/unmapped_bam; \
		$$(HISAT) $$(HISAT_OPTS) \
		--rg-id $2 --rg \"SM:$1\" \
		--rg PL:$${SEQ_PLATFORM} --rg \"LB:\$1\" \
		-S >($$(SAMTOOLS2) view -bh - > hisat/bam/$2.hisat.bam) \
		--un >($$(SAMTOOLS2) view -bh - > hisat/unmapped_bam/$2.hisat_unmapped.bam) \
		$$(if $$(<<),-1 $$(<) -2 $$(<<),-U $$(<))")
hisat/unmapped_bam/$2.hisat_unmapped.bam : hisat/bam/$2.hisat.bam
endef
$(foreach ss,$(SPLIT_SAMPLES),$(eval $(call align-merged-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss)))))
endif


include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
