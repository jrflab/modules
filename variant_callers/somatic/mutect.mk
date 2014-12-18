##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

VCF_SAMPLES = 0 1
VCF_GEN_IDS = GT AD DP FA

SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer -canon

MUTECT_JAR := $(HOME)/share/usr/lib/java/muTect-1.1.4.jar
MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_OPTS = --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION)
MUTECT = $(JAVA) -Xmx11G -jar $(MUTECT_JAR) --analysis_type MuTect $(MUTECT_OPTS)

MUT_FREQ_REPORT = $(RSCRIPT) $(HOME)/share/scripts/plotSeqLogoFromMutect.R


MUTECT_FILTER_SUFFIX := som_ad_ft
ifdef TARGETS_FILE
TARGET_FILTER ?= true
ifeq ($(TARGET_FILTER),true)
MUTECT_FILTER_SUFFIX := $(MUTECT_FILTER_SUFFIX).target_ft
endif
endif
MUTECT_FILTER_SUFFIX := $(MUTECT_FILTER_SUFFIX).pass.dbsnp.nsfp.eff
EFF_TYPES = silent missense nonsilent_cds nonsilent
MUTECT_VCF_SUFFIXES = mutect.$(MUTECT_FILTER_SUFFIX)
MUTECT_TABLE_SUFFIXES = mutect.$(MUTECT_FILTER_SUFFIX).tab $(foreach eff,$(EFF_TYPES),mutect.$(MUTECT_FILTER_SUFFIX).tab.$(eff).novel mutect.$(MUTECT_FILTER_SUFFIX).tab.$(eff))

#VCFS = $(foreach suff,$(VCF_SUFFIXES),$(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).$(suff).vcf))
VCFS = $(foreach suff,$(MUTECT_VCF_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(suff).vcf))

# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect-tumor-normal-chr
mutect/chr_vcf/$1_$2.$3.mutect%vcf mutect/chr_tables/$1_$2.$3.mutect%txt : bam/$1%bam bam/$2%bam
	$$(MKDIR) mutect/chr_tables mutect/chr_vcf; $$(call LSCRIPT_MEM,12G,16G,"$$(MUTECT) --enable_extended_output --intervals $3 --reference_sequence $$(REF_FASTA) --cosmic $$(COSMIC) --dbsnp $$(DBSNP1PC) --input_file:tumor $$< --input_file:normal $$(word 2,$$^) -vcf mutect/chr_vcf/$1_$2.$3.mutect.vcf --out mutect/chr_tables/$1_$2.$3.mutect.txt")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach i,$(SETS_SEQ), \
		$(foreach tumor,$(call get_tumors,$(set.$i)), \
			$(eval $(call mutect-tumor-normal-chr,$(tumor),$(call get_normal,$(set.$i)),$(chr))))))

# merge variant tables 
define ext-mutect-tumor-normal
mutect/tables/$1.mutect.txt : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_tables/$1.$$(chr).mutect.txt)
	$$(INIT) head -2 $$< > $$@; for table in $$^; do sed '1,2d' $$$$table >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ext-mutect-tumor-normal,$(pair))))

mutect/report/report.timestamp: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_MEM,3G,5G,"$(MUT_FREQ_REPORT) --outDir $(@D) $^ && touch $@")

mutect/lowAFreport/report.timestamp: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_MEM,3G,5G,"$(MUT_FREQ_REPORT) --outDir $(@D) --lowAF $^ && touch $@")

mutect/highAFreport/report.timestamp: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_MEM,3G,5G,"$(MUT_FREQ_REPORT) --outDir $(@D) --highAF $^ && touch $@")

# merge variants 
#$$(INIT) grep '^##' $$< > $$@; echo "##PEDIGREE=<Derived=$1,Original=$2>" >> $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
define mutect-tumor-normal
vcf/$1_$2.mutect.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_vcf/$1_$2.$$(chr).mutect.vcf)
	$$(INIT) grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call mutect-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

include ~/share/modules/vcf_tools/vcftools.mk

