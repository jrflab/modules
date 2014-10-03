# Run VarScan
##### DEFAULTS ######


REF ?= hg19
LOGDIR = log/varscan.$(NOW)

SPLIT_CHR ?= true

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

VARSCAN_JAR = $(JARDIR)/VarScan.v2.3.6.jar
VARSCAN_MEM = $(JAVA) -Xmx$1 -jar $(VARSCAN_JAR)
VARSCAN = $(call VARSCAN_MEM,8G)
FIX_VARSCAN_VCF = $(PERL) ~/share/scripts/fixVarscanVcf.pl


VPATH ?= bam

.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: all vcfs copycalls segments cnv reports

SNP_VCF_EFF_FIELDS += VAF
INDEL_VCF_EFF_FIELDS += VAF

ANN_TYPES = eff # annotated
EFF_TYPES = silent missense nonsilent_cds nonsilent
VARIANT_TYPES = varscan_snps varscan_indels

MIN_MAP_QUAL ?= 1
VARSCAN_MIN_COVERAGE ?= 8
VARSCAN_MIN_READS2 ?= 2
VARSCAN_MIN_AVG_QUAL ?= 15
VARSCAN_MIN_VAR_FREQ ?= 0.1
VARSCAN_MIN_FREQ_FOR_HOM ?= 0.75
VARSCAN_P_VALUE ?= 99e-02
VARSCAN_STRAND_FILTER ?= 1

override VARSCAN_OPTS = --min-coverage $(VARSCAN_MIN_COVERAGE) \
	--min-reads2 $(VARSCAN_MIN_READS2) \
	--min-avg-qual $(VARSCAN_MIN_AVG_QUAL) \
	--min-var-freq $(VARSCAN_MIN_VAR_FREQ) \
	--min-freq-for-hom $(VARSCAN_MIN_FREQ_FOR_HOM) \
	--p-value $(VARSCAN_P_VALUE) \
	--strand-filter $(VARSCAN_STRAND_FILTER) 

FILTER_SUFFIX := dp_ft.fp_ft.dgd_ft.encode_ft
ifdef TARGETS_FILE
FILTER_SUFFIX := $(FILTER_SUFFIX).target_ft
endif
ANN_SUFFIX := pass.dbsnp.eff

VCF_SUFFIX.varscan_snps := $(FILTER_SUFFIX).$(ANN_SUFFIX).nsfp.chasm.fathmm.transfic

VCF_SUFFIX.varscan_indels := $(FILTER_SUFFIX)
ifeq ($(HRUN),true)
HRUN_FILTER ?= 1
VCF_SUFFIX.varscan_indels := $(VCF_SUFFIX.varscan_indels).hrun.hrun_ft
endif
VCF_SUFFIX.varscan_indels := $(VCF_SUFFIX.varscan_indels).$(ANN_SUFFIX)

VCF_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(type).$(VCF_SUFFIX.$(type)))
TABLE_SUFFIXES = $(foreach suff,$(VCF_SUFFIXES),$(foreach eff,$(EFF_TYPES),$(suff).tab.$(eff).novel))

VCFS = $(foreach sample,$(SAMPLES),$(foreach suff,$(VCF_SUFFIXES),vcf/$(sample).$(suff).vcf))
TABLES = $(foreach sample,$(SAMPLES),$(foreach suff,$(TABLE_SUFFIXES),tables/$(sample).$(suff).txt))
TABLES += $(foreach suff,$(TABLE_SUFFIXES),alltables/all.$(suff).txt) alltables/all.varscan_snps.dp_ft.target_ft.pass.dbsnp.eff.nsfp.chasm.fathmm.transfic.tab.txt

all : vcfs tables cnv
variants : vcfs tables
cnv : copycalls segments
vcfs : $(VCFS)
tables : $(TABLES)
reports : $(foreach type,$(VARIANT_TYPES),reports/$(type).$(FILTER_SUFFIX).grp)


ifeq ($(SPLIT_CHR),true)
define varscan-chr-type
varscan/chr_vcf/%.$1.$2.vcf : bam/%.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) mpileup2$2 \
	<($$(SAMTOOLS) mpileup -r $1 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--output-vcf $$(VARSCAN_OPTS) --vcf-sample-list $$* | $$(FIX_VARSCAN_VCF) -s $$* > $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach type,snp indel,$(eval $(call varscan-chr-type,$(chr),$(type)))))

define merge-varscan-vcfs
varscan/vcf/$1.%.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1.$$(chr).%.vcf)
	$$(INIT) grep '^##' $$< > $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-varscan-vcfs,$(sample))))

else # no splitting by chr

define varscan-type
varscan/vcf/%.$1.vcf : bam/%.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) mpileup2$1 \
	<($$(SAMTOOLS) mpileup -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--output-vcf --vcf-sample-list $$* $$(VARSCAN_OPTS) | $$(FIX_VARSCAN_VCF) -s $$* > $$@")
endef
$(foreach type,snp indel, $(eval $(call varscan-type,$(chr),$(type))))
endif

vcf/%.varscan_indels.vcf : varscan/vcf/%.indel.vcf
	$(INIT) ln $< $@

vcf/%.varscan_snps.vcf : varscan/vcf/%.snp.vcf
	$(INIT) ln $< $@

include ~/share/modules/variant_callers/gatk.mk
