# This module is for the Tophat aligner
# input: $(SAMPLES) 
# Options: 
# 	PHRED64 = true/false
# 	NO_NOVEL_SPLICING = true/false
# 	NUM_CORES = 4
#	INNER_MATE_DIST = 200
include modules/Makefile.inc

NO_REALN = true
NO_RECAL = true
NO_FILTER = true
DUP_TYPE = none
include modules/aligners/align.inc

LOGDIR = log/tophat.$(NOW)

NUM_CORES ?= 4
NO_NOVEL_SPLICING ?= false 

TOPHAT = $(HOME)/usr/bin/tophat2
TOPHAT_OPTS = --keep-fasta-order --no-sort-bam -G $(GENES_GTF) -p ${NUM_CORES} --tmp-dir $(TMPDIR)/$*

SEQ_PLATFORM ?= illumina

ifeq ($(PHRED64),true)
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
tophat_bams : $(addsuffix .md5,$(BAMS)) $(addsuffix .bai,$(BAMS))
#tophat_unmapped_bams : $(foreach sample,$(SAMPLES),tophat/unmapped_bam/$(sample).unmapped.bam)
#tophat_junctions: $(foreach sample,$(SAMPLES),tophat/junctions/$(sample)_junctions.bed)

bam/%.bam.md5 : tophat/bam/%.tophat.$(BAM_SUFFIX).md5
	$(INIT) ln -f $(<:.md5=) $(@:.md5=) && $(MD5)

tophat/bam/%.tophat.sorted.bam.md5 : tophat/%/accepted_hits.sorted.bam.md5 tophat/%/unmapped.sorted.bam.md5 tophat/%/accepted_hits.sorted.bam.bai tophat/%/unmapped.sorted.bam.bai
	$(call LSCRIPT_MEM,7G,7G,"$(CHECK_MD5) $(SAMTOOLS) merge -f $(@M) $(<M) $(<<M) && $(MD5)")

tophat/%/accepted_hits.bam.md5 tophat/%/unmapped.bam.md5 : fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5
	$(call LSCRIPT_NAMED_PARALLEL_MEM,$*_tophat,4,6G,10G,"$(CHECK_MD5) LBID=`echo \"$*\" | sed 's/_[0-9]\+//'`; \
		$(TOPHAT) $(TOPHAT_OPTS) \
		--rg-id $* --rg-library \"\$${LBID}\" \
		--rg-platform \"${SEQ_PLATFORM}\" --rg-sample \"\$${LBID}\" \
		-o tophat/$* $(BOWTIE_REF) $(<M) $(<<M) && \
		md5sum tophat/$*/accepted_hits.bam > tophat/$*/accepted_hits.bam.md5 && \
		md5sum tophat/$*/unmapped.bam > tophat/$*/unmapped.bam.md5")

ifdef SPLIT_SAMPLES
define merged-bam
ifeq ($(shell echo "$(words $2) > 1" | bc),1)
tophat/bam/$1.header.sam : $$(foreach split,$2,tophat/bam/$$(split).tophat.sorted.bam.md5)
	$$(INIT) $$(SAMTOOLS) view -H $$(<M) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp

tophat/bam/$1.tophat.sorted.bam.md5 : tophat/bam/$1.header.sam $$(foreach split,$2,tophat/bam/$$(split).tophat.sorted.bam.md5)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$(@M) $$(filter %.bam,$$(^M)) && $$(MD5) && $$(RM) $$(^M) $$^")
endif
ifeq ($(shell echo "$(words $2) == 1" | bc),1)
tophat/bam/$1.tophat.sorted.bam.md5 : tophat/bam/$2.tophat.sorted.bam.md5
	$$(INIT) mv $$(<M) $$(@M) && $$(MD5)
endif
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))
endif


include modules/fastq_tools/fastq.mk
include modules/bam_tools/processBam.mk
