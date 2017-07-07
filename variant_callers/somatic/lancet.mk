# lancet variant detection
# vim: set ft=make :

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/lancet.$(NOW)

LANCET_TARGET_ONLY ?= false
NUM_LANCET_CHUNKS = 50
LANCET_CHUNKS = $(shell seq -f "%03g" $(NUM_LANCET_CHUNKS))

LANCET_MIN_COV ?= 5
LANCET_DIR = $(HOME)/share/usr/lancet
LANCET = export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LANCET_DIR)/bamtools-2.3.0/lib/; \
	$(LANCET_DIR)/lancet
LANCET_OPTS = --ref $(REF_FASTA) --cov-thr $(LANCET_MIN_COV)


..DUMMY := $(shell mkdir -p version; echo "$(LANCET) $(LANCET_OPTS) > version/lancet.txt")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all lancet_vcfs # lancet_mafs

lancet : lancet_vcfs #lancet_mafs


lancet_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).lancet_indels.vcf)
#lancet_mafs : $(foreach pair,$(SAMPLE_PAIR),maf/$(pair).lancet_indels.maf)

ifeq ($(LANCET_TARGET_ONLY),true)
lancet/interval_chunk/chunk.timestamp : $(TARGETS_FILE)
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_LANCET_CHUNKS) $<  && touch $@
else
lancet_genome_sliding_window.bed : $(REF_DICT)
	$(INIT) bedtools makewindows -g <(sed 's/.*SN:\([^\s]\+\)\sLN:\([0-9]\+\).*/\1\t\2/' $<) -w 10000 -s 9900 > $@
lancet/interval_chunk/chunk.timestamp : lancet_genome_sliding_window.bed
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_LANCET_CHUNKS) $<  && touch $@
endif


define lancet-interval-chunk
lancet/interval_chunk/chunk$1.bed : lancet/interval_chunk/chunk.timestamp
endef
$(foreach i,$(LANCET_CHUNKS),$(eval $(call lancet-interval-chunk,$i)))

define lancet-chunk-tumor-normal
lancet/chunk_vcf/$2_$3.lancet.$1.vcf : lancet/interval_chunk/chunk$1.bed bam/$2.bam bam/$3.bam
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(LANCET) --num-threads 4 --tumor $$(<<) --normal $$(<<<) $$(LANCET_OPTS) --bed $$(<) > $$@")
endef
$(foreach chunk,$(LANCET_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lancet-chunk-tumor-normal,$(chunk),$(tumor.$(pair)),$(normal.$(pair))))))

INDEL_FILTER_VCF = python modules/vcf_tools/indel_filter_vcf.py
SNP_FILTER_VCF = python modules/vcf_tools/snp_filter_vcf.py
define lancet-tumor-normal
vcf/$1_$2.lancet_snps.vcf : $$(foreach chunk,$$(LANCET_CHUNKS),lancet/chunk_vcf/$1_$2.lancet.$$(chunk).vcf)
	$$(call LSCRIPT_MEM,4G,8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | \
		$$(VCF_SORT) $$(REF_DICT) -) | $$(SNP_FILTER_VCF) > $$@")
vcf/$1_$2.lancet_indels.vcf : $$(foreach chunk,$$(LANCET_CHUNKS),lancet/chunk_vcf/$1_$2.lancet.$$(chunk).vcf)
	$$(call LSCRIPT_MEM,4G,8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(INDEL_FILTER_VCF) > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lancet-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk
include modules/bam_tools/processBam.mk
