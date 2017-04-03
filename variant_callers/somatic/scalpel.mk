# scalpel variant detection
# vim: set ft=make :

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/scalpel.$(NOW)

NUM_SCALPEL_CHUNKS = 2
SCALPEL_CHUNKS = $(shell seq -f "%03g" $(NUM_SCALPEL_CHUNKS))

SCALPEL_MIN_COV ?= 5
SCALPEL_DIR = $(HOME)/share/usr/scalpel-0.5.3
SCALPEL_DISCOVERY = export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(SCALPEL_DIR)/bamtools-2.3.0/lib/; $(PERL) $(SCALPEL_DIR)/scalpel-discovery
SCALPEL_OPTS = --ref $(REF_FASTA) --format vcf --covthr $(SCALPEL_MIN_COV)


SCALPEL2VCF = $(PERL) modules/variant_callers/somatic/scalpelToVcf.pl

..DUMMY := $(shell mkdir -p version; echo "$(SCALPEL) $(SCALPEL_OPTS) > version/scalpel.txt")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all scalpel_vcfs # scalpel_mafs

scalpel : scalpel_vcfs #scalpel_mafs


scalpel_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).scalpel_indels.vcf)
#scalpel_mafs : $(foreach pair,$(SAMPLE_PAIR),maf/$(pair).scalpel_indels.maf)

ifdef TARGETS_FILE
scalpel/interval_chunk/chunk.timestamp : $(TARGETS_FILE)
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_SCALPEL_CHUNKS) $<  && touch $@
genome_sliding_window.bed : $(REF_DICT)
	$(INIT) bedtools makewindows -g <('s/.*SN:\([^\s]\+\)\sLN:\([0-9]\+\).*/\1\t\2/' $<) -w 10000 -s 9900 > $@
else
scalpel/interval_chunk/chunk.timestamp : genome_sliding_window.bed
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk $(NUM_SCALPEL_CHUNKS) $<  && touch $@
endif


define scalpel-interval-chunk
scalpel/interval_chunk/chunk$1.bed : scalpel/interval_chunk/chunk.timestamp
endef
$(foreach i,$(SCALPEL_CHUNKS),$(eval $(call scalpel-interval-chunk,$i)))

define scalpel-chunk-tumor-normal
scalpel/$2_$3/$1/scalpel.timestamp : scalpel/interval_chunk/chunk$1.bed bam/$2.bam bam/$3.bam
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$2_$3_$1_scalpel,4,2G,3G,"$$(SCALPEL_DISCOVERY) --somatic --numprocs 4 --tumor $$(<<) --normal $$(<<<) $$(SCALPEL_OPTS) --bed $$(<) --dir $$(@D) && touch $$@")
scalpel/$2_$3/$1/somatic.indel.vcf: scalpel/$2_$3/$1/scalpel.timestamp 
	if [ ! -e $$@ ]; then cp $$(@D)/main/$$(@F) $$@; fi
endef
$(foreach chunk,$(SCALPEL_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call scalpel-chunk-tumor-normal,$(chunk),$(tumor.$(pair)),$(normal.$(pair))))))

define scalpel-tumor-normal
vcf/$1_$2.scalpel_indels.vcf : $$(foreach chunk,$$(SCALPEL_CHUNKS),scalpel/$1_$2/$$(chunk)/somatic.indel.vcf)
	$$(call LSCRIPT_MEM,4G,8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk
include modules/bam_tools/processBam.mk
