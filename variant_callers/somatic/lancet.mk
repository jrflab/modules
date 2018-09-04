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

LANCET_THREADS = 4
LANCET_RUN_OPTS = -c -n $(LANCET_THREADS) -s 2G -m 3G -w 36:00:00
LANCET_OPTS = --ref $(REF_FASTA) --cov-thr $(LANCET_MIN_COV) --num-threads $(LANCET_THREADS)

LANCET_SOURCE_ANN_VCF = python modules/vcf_tools/annotate_source_vcf.py --source lancet
TUMOR_VARIANT_READ_FILTER_VCF = python modules/vcf_tools/tumor_variant_read_filter_vcf.py --pass_only

..DUMMY := $(shell mkdir -p version; echo "$(LANCET) $(LANCET_OPTS) > version/lancet.txt")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all lancet_vcfs

lancet : lancet_vcfs


lancet_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).lancet_indels.vcf)

ifeq ($(LANCET_TARGET_ONLY),true)
lancet/interval_chunk/chunk.timestamp : $(TARGETS_FILE)
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_LANCET_CHUNKS) $<  && touch $@

define lancet-interval-chunk
lancet/interval_chunk/chunk$1.bed : lancet/interval_chunk/chunk.timestamp
endef
$(foreach i,$(LANCET_CHUNKS),$(eval $(call lancet-interval-chunk,$i)))

define lancet-chunk-tumor-normal
lancet/chunk_vcf/$2_$3.lancet.$1.vcf : lancet/interval_chunk/chunk$1.bed bam/$2.bam bam/$3.bam
	$$(call RUN,$$(LANCET_RUN_OPTS),"$$(LANCET) --tumor $$(<<) --normal $$(<<<) $$(LANCET_OPTS) \
		--bed $$(<) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach chunk,$(LANCET_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lancet-chunk-tumor-normal,$(chunk),$(tumor.$(pair)),$(normal.$(pair))))))

define lancet-tumor-normal
lancet/vcf/$1_$2.lancet_snps.vcf : $$(foreach chunk,$$(LANCET_CHUNKS),lancet/chunk_vcf/$1_$2.lancet.$$(chunk).vcf)
	$$(call RUN,-s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(SNP_FILTER_VCF) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
lancet/vcf/$1_$2.lancet_indels.vcf : $$(foreach chunk,$$(LANCET_CHUNKS),lancet/chunk_vcf/$1_$2.lancet.$$(chunk).vcf)
	$$(call RUN,-s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(INDEL_FILTER_VCF) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lancet-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

else
define lancet-chr-tumor-normal
lancet/chr_vcf/$2_$3.lancet.$1.vcf : bam/$2.bam bam/$3.bam
	$$(call RUN,$$(LANCET_RUN_OPTS),"$$(LANCET) --tumor $$(<) --normal $$(<<) \
		$$(LANCET_OPTS) --reg $1 > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lancet-chr-tumor-normal,$(chr),$(tumor.$(pair)),$(normal.$(pair))))))

define lancet-tumor-normal
lancet/vcf/$1_$2.lancet_snps.vcf : $$(foreach chr,$$(CHROMOSOMES),lancet/chr_vcf/$1_$2.lancet.$$(chr).vcf)
	$$(call RUN,-s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(SNP_FILTER_VCF) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
lancet/vcf/$1_$2.lancet_indels.vcf : $$(foreach chr,$$(CHROMOSOMES),lancet/chr_vcf/$1_$2.lancet.$$(chr).vcf)
	$$(call RUN,-s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(INDEL_FILTER_VCF) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call lancet-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

define filter_lancet-tumor-normal
vcf/$1_$2.lancet_%.vcf lancet/vcf/$1_$2.lancet_%.vcf.tmp : lancet/vcf/$1_$2.lancet_%.vcf
	$$(call RUN,-s 1G -m 2G,"$$(RSCRIPT) modules/scripts/swapvcf.R --file $$< --tumor $1 --normal $2 && \
							 cp $$<.tmp $$< && \
							 $$(LANCET_SOURCE_ANN_VCF) < $$< | \
							 $$(TUMOR_VARIANT_READ_FILTER_VCF) -t $1 -n $2 > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call filter_lancet-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk
include modules/bam_tools/processBam.mk
