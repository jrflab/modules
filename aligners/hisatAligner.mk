include modules/Makefile.inc

BAM_NO_REALN = true
BAM_NO_RECAL = true
BAM_NO_FILTER = true
BAM_DUP_TYPE = none

ALIGNER := hisat

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
	$(call RUN,-N $*_hisat -n 8 -s 1G -m 1.5G,"mkdir -p hisat/bam hisat/unmapped_bam; \
		LBID=`echo \"$*\" | sed 's/_.\+//'`; \
		$(HISAT) $(HISAT_OPTS) \
		--rg-id $* --rg \"SM:\$${LBID}\" \
		--rg PL:${SEQ_PLATFORM} --rg \"LB:\$${LBID}\" \
		-S >($(SAMTOOLS2) view -bh - > hisat/bam/$*.hisat.bam) \
		--un >($(SAMTOOLS2) view -bh - > hisat/unmapped_bam/$*.hisat_unmapped.bam) \
		-1 $(<) -2 $(<<)")

# usage $(eval $(call align-merged-split-fastq,merged-sample,split,fastq(s)))
define align-merged-split-fastq
hisat/bam/$2.hisat.bam : $3
	$$(call RUN,-N $2_hisat -n 8 -s 1G -m 1.5G,"mkdir -p hisat/bam hisat/unmapped_bam; \
		$$(HISAT) $$(HISAT_OPTS) \
		--rg-id $2 --rg \"SM:$1\" \
		--rg PL:$${SEQ_PLATFORM} --rg \"LB:\$1\" \
		-S >($$(SAMTOOLS2) view -bh - > hisat/bam/$2.hisat.bam) \
		--un >($$(SAMTOOLS2) view -bh - > hisat/unmapped_bam/$2.hisat_unmapped.bam) \
		$$(if $$(<<),-1 $$(<) -2 $$(<<),-U $$(<))")
hisat/unmapped_bam/$2.hisat_unmapped.bam : hisat/bam/$2.hisat.bam
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-merged-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))


include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
include modules/aligners/align.mk
