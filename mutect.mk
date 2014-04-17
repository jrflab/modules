# Run mutect on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/mutect.$(NOW)
VCF_SAMPLES = 0 1
VCF_GEN_IDS = GT AD DP FA

SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer -canon

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

MUTECT_JAR := $(HOME)/share/usr/lib/java/muTect-1.1.4.jar
MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_OPTS = --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION)
MUTECT = $(JAVA) -Xmx11G -jar $(MUTECT_JAR) --analysis_type MuTect $(MUTECT_OPTS)

MUT_FREQ_REPORT = $(RSCRIPT) $(HOME)/share/scripts/plotSeqLogoFromMutect.R

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all mutect_vcfs mutect_tables ext_output mut_report

all : mutect_vcfs mutect_tables ext_output mut_report


FILTER_SUFFIX := som_ad_ft.pass.dbsnp.nsfp.eff.chasm.fathmm.transfic
EFF_TYPES = silent missense nonsilent_cds nonsilent
VCF_SUFFIXES = mutect.$(FILTER_SUFFIX)
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),mutect.$(FILTER_SUFFIX).tab.$(eff).novel)

#VCFS = $(foreach suff,$(VCF_SUFFIXES),$(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).$(suff).vcf))
VCFS = $(foreach suff,$(VCF_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(suff).vcf))
mutect_vcfs : $(VCFS) $(addsuffix .idx,$(VCFS))
mutect_tables : $(foreach suff,$(TABLE_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).$(suff).txt)) \
	$(foreach suff,$(TABLE_SUFFIXES),alltables/allTN.$(suff).txt)
ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
mut_report : mutect/report/report.timestamp mutect/lowAFreport/report.timestamp mutect/highAFreport/report.timestamp

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
define mutect-tumor-normal
vcf/$1_$2.mutect.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_vcf/$1_$2.$$(chr).mutect.vcf)
	$$(INIT) grep '^##' $$< > $$@; echo "##PEDIGREE=<Derived=$1,Original=$2>" >> $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call mutect-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))


include ~/share/modules/vcftools.mk

