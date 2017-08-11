# Run mutect on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_FILTERS = DuplicateRead FailsVendorQualityCheck NotPrimaryAlignment BadMate MappingQualityUnavailable UnmappedRead BadCigar
MUTECT_OPTS ?= --enable_extended_output --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) \
			   --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION) -R $(REF_FASTA) \
			   --dbsnp $(DBSNP) $(foreach ft,$(MUTECT_FILTERS),-rf $(ft))
MUTECT = $(JAVA7) -Xmx11G -jar $(MUTECT_JAR) --analysis_type MuTect $(MUTECT_OPTS)

MUTECT_USE_CONTEST ?= false

MUTECT_SPLIT_CHR ?= true

MUT_FREQ_REPORT = modules/variant_callers/somatic/mutectReport.Rmd

MUTECT_SOURCE_ANN_VCF = python modules/vcf_tools/annotate_source_vcf.py --source mutect

..DUMMY := $(shell mkdir -p version; echo "$(MUTECT) &> version/mutect.txt")

LOGDIR ?= log/mutect.$(NOW)

PHONY += mutect mutect_vcfs mutect_mafs mut_report
mutect : mutect_vcfs # mutect_mafs
mutect_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).mutect.vcf)
mut_report : mutect/report/index.html mutect/lowAFreport/index.html mutect/highAFreport/index.html

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)

SPLIT_BED = python modules/scripts/split_bed.py
NUM_MUTECT_CHUNKS = 100
MUTECT_CHUNKS = $(shell seq -w 1 $(NUM_MUTECT_CHUNKS))

ifdef TARGETS_FILE
mutect/interval_chunk/chunk.timestamp : $(TARGETS_FILE)
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk --num_chunks $(NUM_MUTECT_CHUNKS) $<  && touch $@
else
genome_sliding_window.bed : $(REF_DICT)
	$(INIT) bedtools makewindows -g <('s/.*SN:\([^\s]\+\)\sLN:\([0-9]\+\).*/\1\t\2/' $<) -w 10000 -s 9900 > $@
mutect/interval_chunk/chunk.timestamp : genome_sliding_window.bed
	$(INIT) $(SPLIT_BED) --out_prefix $(@D)/chunk $(NUM_MUTECT_CHUNKS) $<  && touch $@
endif

define interval-chunk
mutect/interval_chunk/chunk$1.bed : mutect/interval_chunk/chunk.timestamp
endef
$(foreach i,$(MUTECT_CHUNKS),$(eval $(call interval-chunk,$i)))

# run mutect on each chunk
#$(call mutect-tumor-normal-chunk,tumor,normal,chunk)
define mutect-tumor-normal-chunk
mutect/chunk_vcf/$1_$2.chunk$3.mutect.vcf : mutect/interval_chunk/chunk$3.bed bam/$1.bam bam/$2.bam \
	$$(if $$(findstring true,$$(MUTECT_USE_CONTEST)),contest/$1_$2.contest.txt) bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,-c -s 12G -m 15G,"$$(MKDIR) mutect/chunk_tables mutect/cov; \
		$$(MUTECT) --tumor_sample_name $1 --normal_sample_name $2 \
		$$(if $$(findstring true,$$(MUTECT_USE_CONTEST)),--fraction_contamination `csvcut -t -c contamination $$(<<<<) | sed 1d`) \
		--intervals $$< \
		-I:tumor $$(<<) -I:normal $$(<<<) --out mutect/chunk_tables/$1_$2.chunk$3.mutect.txt -vcf $$@ \
		--coverage_file mutect/cov/$1_$2.chunk$3.cov.wig")
endef
$(foreach chunk,$(MUTECT_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal-chunk,$(tumor.$(pair)),$(normal.$(pair)),$(chunk)))))

# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect-tumor-normal-chr
mutect/chr_vcf/$1_$2.$3.mutect.vcf : bam/$1.bam bam/$2.bam \
	$$(if $$(findstring true,$$(MUTECT_USE_CONTEST)),contest/$1_$2.contest.txt) bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,-c -s 12G -m 15G,"$$(MKDIR) mutect/chr_tables mutect/cov; \
		$$(MUTECT) --tumor_sample_name $1 --normal_sample_name $2 \
		$$(if $$(findstring true,$$(MUTECT_USE_CONTEST)),--fraction_contamination `csvcut -t -c contamination $$(<<<) | sed 1d`) \
		--intervals $3 -I:tumor $$(<) -I:normal $$(<<) --out mutect/chr_tables/$1_$2.$3.mutect.txt  \
		-vcf $$@ --coverage_file mutect/cov/$1_$2.$3.cov.wig")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))


# merge variant tables 
ifeq ($(MUTECT_SPLIT_CHR),true)
define ext-mutect-tumor-normal
mutect/tables/$1.mutect.txt : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_tables/$1.$$(chr).mutect.txt)
	$$(INIT) head -2 $$< > $$@; for table in $$^; do sed '1,2d' $$$$table >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ext-mutect-tumor-normal,$(pair))))

define mutect-tumor-normal
vcf/$1_$2.mutect.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_vcf/$1_$2.$$(chr).mutect.vcf)
	$$(call RUN,-c -s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | grep PASS | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(MUTECT_SOURCE_ANN_VCF) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
else
define ext-mutect-tumor-normal
mutect/tables/$1.mutect.txt : $$(foreach chunk,$$(MUTECT_CHUNKS),mutect/chunk_tables/$1.chunk$$(chunk).mutect.txt)
	$$(INIT) head -2 $$< > $$@; for table in $$^; do sed '1,2d' $$$$table >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ext-mutect-tumor-normal,$(pair))))

define mutect-tumor-normal
vcf/$1_$2.mutect.vcf : $$(foreach chunk,$$(MUTECT_CHUNKS),mutect/chunk_vcf/$1_$2.chunk$$(chunk).mutect.vcf)
	$$(call RUN,-c -s 4G -m 8G,"(grep '^#' $$<; cat $$^ | grep -v '^#' | grep PASS | $$(VCF_SORT) $$(REF_DICT) -) | \
		$$(MUTECT_SOURCE_ANN_VCF) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

mutect/report/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call RUN,-N mutect_report -s 6G -m 35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) $^")

mutect/lowAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call RUN,-N mutect_lowaf_report -s 6G -m 35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --lowAF $^")

mutect/highAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call RUN,-N mutect_highaf_report -s 6G -m 35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --highAF $^")

include modules/vcf_tools/vcftools.mk
include modules/contamination/contest.mk
