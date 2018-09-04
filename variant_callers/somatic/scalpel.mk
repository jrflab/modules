# scalpel variant detection
# vim: set ft=make :

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/scalpel.$(NOW)

SCALPEL_TARGET_ONLY ?= true

NUM_SCALPEL_CHUNKS = 30
SCALPEL_CHUNKS = $(shell seq -f "%03g" $(NUM_SCALPEL_CHUNKS))

SCALPEL_MIN_COV ?= 5
SCALPEL_DIR = $(HOME)/share/usr/scalpel-0.5.3
SCALPEL_DISCOVERY = export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(SCALPEL_DIR)/bamtools-2.3.0/lib/; $(PERL) $(SCALPEL_DIR)/scalpel-discovery
SCALPEL_OPTS = --ref $(REF_FASTA) --format vcf --covthr $(SCALPEL_MIN_COV)

SCALPEL2VCF = $(PERL) modules/variant_callers/somatic/scalpelToVcf.pl

SCALPEL_SOURCE_ANN_VCF = python modules/vcf_tools/annotate_source_vcf.py --source scalpel
TUMOR_VARIANT_READ_FILTER_VCF = python modules/vcf_tools/tumor_variant_read_filter_vcf.py --pass_only

..DUMMY := $(shell mkdir -p version; echo "$(SCALPEL) $(SCALPEL_OPTS) > version/scalpel.txt")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all scalpel_vcfs # scalpel_mafs

scalpel : scalpel_vcfs #scalpel_mafs

scalpel_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).scalpel_indels.vcf)
#scalpel_mafs : $(foreach pair,$(SAMPLE_PAIR),maf/$(pair).scalpel_indels.maf)

ifeq ($(SCALPEL_TARGET_ONLY),true)
scalpel/interval_chunk/chunk.timestamp : $(TARGETS_FILE)
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_SCALPEL_CHUNKS) $<  && touch $@
#else
#scalpel_genome_sliding_window.bed : $(REF_DICT)
	#$(INIT) bedtools makewindows -g <(sed 's/.*SN:\([^\s]\+\)\sLN:\([0-9]\+\).*/\1\t\2/' $<) -w 10000 -s 9900 > $@
#scalpel/interval_chunk/chunk.timestamp : scalpel_genome_sliding_window.bed
	#$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_SCALPEL_CHUNKS) $<  && touch $@
#endif


define scalpel-interval-chunk
scalpel/interval_chunk/chunk$1.bed : scalpel/interval_chunk/chunk.timestamp
endef
$(foreach i,$(SCALPEL_CHUNKS),$(eval $(call scalpel-interval-chunk,$i)))

define scalpel-chunk-tumor-normal
scalpel/$2_$3/$1/main/somatic.indel.vcf : scalpel/interval_chunk/chunk$1.bed bam/$2.bam bam/$3.bam
	$$(call RUN,-N $2_$3_$1_scalpel -c -n 4 -s 2G -m 3G,"$$(SCALPEL_DISCOVERY) --somatic --numprocs 4 --tumor $$(<<) --normal $$(<<<) $$(SCALPEL_OPTS) --bed $$(<) --dir $$(@D)/..")
endef
$(foreach chunk,$(SCALPEL_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call scalpel-chunk-tumor-normal,$(chunk),$(tumor.$(pair)),$(normal.$(pair))))))

define scalpel-tumor-normal
vcf/$1_$2.scalpel_indels.vcf : $$(foreach chunk,$$(SCALPEL_CHUNKS),scalpel/$1_$2/$$(chunk)/main/somatic.indel.vcf)
	$$(call RUN,-c -s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(SCALPEL_SOURCE_ANN_VCF) | $$(TUMOR_VARIANT_READ_FILTER_VCF) -t $1 -n $2 > $$@.tmp && \
		$$(call VERIFY_VCF,$$@.tmp,$$@) && \
		$$(RSCRIPT) modules/scripts/swapvcf.R --file $$@ --tumor $1 --normal $2 && \
		$$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
else

define scalpel-chr-tumor-normal
scalpel/$2_$3/$1/main/somatic.indel.vcf : bam/$2.bam bam/$3.bam
	$$(call RUN,-N $2_$3_$1_scalpel -c -n 4 -s 2G -m 3G,"$$(SCALPEL_DISCOVERY) --somatic --numprocs 4 \
		--tumor $$(<) --normal $$(<<) $$(SCALPEL_OPTS) --bed $1 --dir $$(@D)/..")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call scalpel-chr-tumor-normal,$(chr),$(tumor.$(pair)),$(normal.$(pair))))))

define scalpel-tumor-normal
vcf/$1_$2.scalpel_indels.vcf : $$(foreach chr,$$(CHROMOSOMES),scalpel/$1_$2/$$(chr)/main/somatic.indel.vcf)
	$$(call RUN,-c -s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(SCALPEL_SOURCE_ANN_VCF) | $$(TUMOR_VARIANT_READ_FILTER_VCF) -t $1 -n $2 > $$@.tmp && \
		$$(call VERIFY_VCF,$$@.tmp,$$@) && \
		$$(RSCRIPT) modules/scripts/swapvcf.R --file $$@ --tumor $1 --normal $2 && \
		$$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

endif

include modules/vcf_tools/vcftools.mk
include modules/bam_tools/processBam.mk
