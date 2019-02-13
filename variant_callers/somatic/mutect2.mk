#### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

MUTECT2_MAX_ALT_IN_NORMAL ?= 500
MUTECT2_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT2_FILTERS = DuplicateRead FailsVendorQualityCheck NotPrimaryAlignment BadMate MappingQualityUnavailable UnmappedRead BadCigar
MUTECT2_OPTS ?= --max_alt_alleles_in_normal_count $(MUTECT2_MAX_ALT_IN_NORMAL) --max_alt_allele_in_normal_fraction $(MUTECT2_MAX_ALT_IN_NORMAL_FRACTION) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach ft,$(MUTECT2_FILTERS),-rf $(ft))
MUTECT2 = $(call GATK_MEM2,10G) -T MuTect2 $(MUTECT2_OPTS)

SPLIT_SNPS_INDELS_VCF = python modules/vcf_tools/split_snps_indels_vcf.py

LOGDIR ?= log/mutect2.$(NOW)

PHONY += mutect2 mutect2_vcfs mutect2_tables

mutect2 : mutect2_vcfs mutect2_tables

..DUMMY := $(shell mkdir -p version; echo "$(MUTECT2) &> version/mutect2.txt")

MUTECT2_VARIANT_TYPES = mutect_snps mutect_indels
mutect2_vcfs : $(foreach type,$(MUTECT2_VARIANT_TYPES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(type).vcf))

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)

SPLIT_BED = python modules/scripts/split_bed.py
NUM_MUTECT2_CHUNKS = 50
MUTECT2_CHUNKS = $(shell seq -f "%03g" $(NUM_MUTECT2_CHUNKS))

ifdef TARGETS_FILE
mutect2/interval_chunk/chunk.timestamp : $(TARGETS_FILE)
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_MUTECT2_CHUNKS) $<  && touch $@
else
genome_sliding_window.bed : $(REF_DICT)
	$(INIT) bedtools makewindows -g <('s/.*SN:\([^\s]\+\)\sLN:\([0-9]\+\).*/\1\t\2/' $<) -w 10000 -s 9900 > $@
mutect2/interval_chunk/chunk.timestamp : genome_sliding_window.bed
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk $(NUM_MUTECT2_CHUNKS) $<  && touch $@
endif

define interval-chunk
mutect2/interval_chunk/chunk$1.bed : mutect2/interval_chunk/chunk.timestamp
endef
$(foreach i,$(MUTECT2_CHUNKS),$(eval $(call interval-chunk,$i)))

# run mutect on each chromosome
# $(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect2-tumor-normal-chunk
mutect2/chunk_vcf/$1_$2.chunk$3.mutect_snps_indels.vcf.gz : mutect2/interval_chunk/chunk$3.bed bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,-c -s 12G -m 12G,"$$(MUTECT2) --intervals $$< -I:tumor $$(<<) -I:normal $$(<<<) | bgzip -c > $$@ ")
endef
$(foreach chunk,$(MUTECT2_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect2-tumor-normal-chunk,$(tumor.$(pair)),$(normal.$(pair)),$(chunk)))))

define mutect2-tumor-normal
mutect2/vcf/$1_$2.mutect_snps_indels.vcf : $$(foreach chunk,$$(MUTECT2_CHUNKS),mutect2/chunk_vcf/$1_$2.chunk$$(chunk).mutect_snps_indels.vcf.gz \
	mutect2/chunk_vcf/$1_$2.chunk$$(chunk).mutect_snps_indels.vcf.gz.tbi)
	$$(call RUN,-s 4G -m 8G,"$$(BCFTOOLS) concat -D -a -o $$@ -O v $$(filter %.vcf.gz,$$^)")

vcf/$1_$2.%.vcf : mutect2/vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect2-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

mutect2/vcf/%.mutect_snps.vcf mutect2/vcf/%.mutect_indels.vcf: mutect2/vcf/%.mutect_snps_indels.vcf
	$(call RUN,-s 4G -m 6G,"$(SPLIT_SNPS_INDELS_VCF) -s mutect2/vcf/$*.mutect_snps.vcf -i mutect2/vcf/$*.mutect_indels.vcf $<")

include modules/vcf_tools/vcftools.mk
