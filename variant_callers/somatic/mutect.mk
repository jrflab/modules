# Run mutect on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_FILTERS = DuplicateRead FailsVendorQualityCheck NotPrimaryAlignment BadMate MappingQualityUnavailable UnmappedRead BadCigar
MUTECT_OPTS ?= --enable_extended_output --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) \
			   --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach ft,$(MUTECT_FILTERS),-rf $(ft))
MUTECT = $(JAVA7) -Xmx11G -jar $(MUTECT_JAR) --analysis_type MuTect $(MUTECT_OPTS)

MUTECT_SPLIT_CHR ?= false

MUT_FREQ_REPORT = modules/variant_callers/somatic/mutectReport.Rmd

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
mutect/chunk_vcf/$1_$2.chunk$3.mutect.vcf : mutect/interval_chunk/chunk$3.bed bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_CHECK_MEM,12G,15G,"$$(MKDIR) mutect/chunk_tables; \
		$$(MUTECT) --tumor_sample_name $1 --normal_sample_name $2 --intervals $$< -I:tumor $$(<<) -I:normal $$(<<<) --out mutect/chunk_tables/$1_$2.chunk$3.mutect.txt -vcf $$@")
endef
$(foreach chunk,$(MUTECT_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal-chunk,$(tumor.$(pair)),$(normal.$(pair)),$(chunk)))))

# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect-tumor-normal-chr
mutect/chr_vcf/$1_$2.$3.mutect.vcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_CHECK_MEM,12G,15G,"$$(MKDIR) mutect/chr_tables; \
		$$(MUTECT) --tumor_sample_name $1 --normal_sample_name $2 --intervals $3 -I:tumor $$(<) -I:normal $$(<<) --out mutect/chr_tables/$1_$2.$3.mutect.txt -vcf $$@")
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
	$$(call LSCRIPT_MEM,4G,8G,"grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@")
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
	$$(call LSCRIPT_MEM,4G,8G,"grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

mutect/report/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_NAMED_MEM,mutect_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) $^")

mutect/lowAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_NAMED_MEM,mutect_lowaf_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --lowAF $^")

mutect/highAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_NAMED_MEM,mutect_highaf_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --highAF $^")

include modules/vcf_tools/vcftools.mk
