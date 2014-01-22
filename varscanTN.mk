# Run VarScan on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######


REF ?= hg19
LOGDIR = log/varscan.$(NOW)

SPLIT_CHR ?= true

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

VARSCAN_JAR = $(JARDIR)/VarScan.v2.3.6.jar
VARSCAN_MEM = $(JAVA) -Xmx$1 -jar $(VARSCAN_JAR)
VARSCAN = $(call VARSCAN_MEM,8G)
SEGMENTCNV = $(HOME)/share/scripts/segmentCNV2.R

MIN_MAP_QUAL ?= 1

MIN_VAR_FREQ ?= 0.05

VPATH ?= bam

.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: all vcfs copycalls segments cnv

SNP_VCF_EFF_FIELDS += FREQ
INDEL_VCF_EFF_FIELDS += FREQ

ANN_TYPES = eff # annotated
EFF_TYPES = silent missense nonsilent_cds nonsilent
VARIANT_TYPES = varscan_snps varscan_indels

FILTER_SUFFIX := vdp_ft.freq_ft.dbsnp
FILTER_SUFFIX.varscan_snps := $(FILTER_SUFFIX).nsfp.chasm.fathmm
FILTER_SUFFIX.varscan_indels := $(FILTER_SUFFIX)
VCF_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(foreach ann,$(ANN_TYPES),$(type).$(FILTER_SUFFIX.$(type)).$(ann)))
TABLE_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(foreach ann,$(ANN_TYPES),$(foreach eff,$(EFF_TYPES),$(type).$(FILTER_SUFFIX.$(type)).$(ann).tab.$(eff).pass.novel)))

VCFS = $(foreach pair,$(SAMPLE_PAIRS),$(foreach suff,$(VCF_SUFFIXES),vcf/$(pair).$(suff).vcf))
TABLES = $(foreach pair,$(SAMPLE_PAIRS),$(foreach suff,$(TABLE_SUFFIXES),tables/$(pair).$(suff).txt))
TABLES += $(foreach suff,$(TABLE_SUFFIXES),alltables/allTN.$(suff).txt)

all : vcfs tables cnv
variants : vcfs tables
cnv : copycalls segments
vcfs : $(VCFS)
tables : $(TABLES)
copycalls : $(foreach pair,$(SAMPLE_PAIRS),varscan/copycall/$(pair).copycall)
segments : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).varscan2copynumber.txt)


%.Somatic.vcf : %.vcf
	$(call LSCRIPT_MEM,5G,8G,"$(call VARSCAN_MEM,4G) somaticFilter $< && $(call VARSCAN_MEM,4G) processSomatic $<")

ifeq ($(SPLIT_CHR),true)
#varscan/chr_vcf/$1_$2.$3.indels.vcf varscan/chr_vcf/$1_$2.$3.snps.vcf : $1.bam $2.bam
define varscan-somatic-tumor-normal-chr
varscan/chr_vcf/$1_$2.$3.varscan_timestamp : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) somatic \
	<($$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$(word 2,$$^)) \
	<($$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--min-var-freq $(MIN_VAR_FREQ) --output-vcf --output-indel varscan/chr_vcf/$1_$2.$3.indels.vcf --output-snp varscan/chr_vcf/$1_$2.$3.snps.vcf && touch $$@")

varscan/chr_vcf/$1_$2.$3.indels.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp

varscan/chr_vcf/$1_$2.$3.snps.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp

endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach i,$(SETS_SEQ), \
		$(foreach tumor,$(call get_tumors,$(set.$i)), \
			$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor),$(call get_normal,$(set.$i)),$(chr))))))
#$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor),$(normal_lookup.$(tumor)),$(chr)))))


define varscan-ped-tumor-normal
varscan/vcf/$1_$2.$3.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1_$2.$$(chr).$3.vcf)
	$$(INIT) echo "##PEDIGREE=<Derived=$1,Original=$2>" > $$@ && \
	grep '^#' $$< >> $$@ && \
	cat $$^ | grep -v '^#' |  perl -lane '$$$$F[3] =~ s:/:,:; $$$$F[4] =~ s:/:,:; unless ($$$$F[3] =~ /[+-,]/ || $$$$F[4] =~ /[+-]/) { print join "\t", @F; }' | $$(VCF_SORT) $$(REF_DICT) -  >> $$@
endef
$(foreach type,indels snps, \
	$(foreach i,$(SETS_SEQ),\
		$(foreach tumor,$(call get_tumors,$(set.$i)), \
			$(eval $(call varscan-ped-tumor-normal,$(tumor),$(call get_normal,$(set.$i)),$(type))))))

vcf/%.varscan_indels.vcf : varscan/vcf/%.indels.Somatic.vcf
	$(INIT) ln $< $@

vcf/%.varscan_snps.vcf : varscan/vcf/%.snps.Somatic.vcf
	$(INIT) ln $< $@

else # no splitting by chr

define varscan-somatic-tumor-normal
varscan/vcf/$1_$2.varscan_timestamp : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) somatic \
	<($$(SAMTOOLS) mpileup -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$(word 2,$$^)) \
	<($$(SAMTOOLS) mpileup -q $$(MIN_MAP_QUAL) $$(REF_FASTA) $$<) \
	--min-var-freq $(MIN_VAR_FREQ) --output-vcf --output-indel varscan/vcf/$1_$2.indels.vcf --output-snp varscan/vcf/$1_$2.snps.vcf")

varscan/vcf/$1_$2.indels.vcf : varscan/vcf/$1_$2.varscan_timestamp

varscan/vcf/$1_$2.snps.vcf : varscan/vcf/$1_$2.varscan_timestamp

vcf/$1_$2.varscan_%.vcf : varscan/vcf/$1_$2.%.Somatic.vcf
	$$(INIT) echo "##PEDIGREE=<Derived=$1,Original=$2>" > $$@ && cat $$< >> $$@
endef
#$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call varscan-somatic-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call varscan-somatic-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))
endif

define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.copynumber :  $1.bam $2.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(SAMTOOLS) mpileup -q 1 -f $$(REF_FASTA) $$(word 2,$$^) $$< | awk 'NF == 9 { print }' |  $$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1 &> $$(LOG)")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call varscan-copynum-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

varscan/copycall/%.copycall : varscan/copynum/%.copynumber
	$(call LSCRIPT_MEM,9G,12G,"n=`awk '{ total += \$$7 } END { print total / NR }' $<`; \
	if [ \$$(bc <<< \"\$$n > 0\") -eq 1 ]; then \
		recenter_opt=\"--recenter-up \$$n\"; \
	else \
		n=\$$(bc <<< \"\$$n*-1\"); \
		recenter_opt=\"--recenter-down \$$n\"; \
	fi; \
	$(VARSCAN) copyCaller $< --output-file $@ \$$recenter_opt")

varscan/segment/%.varscan2copynumber.txt : varscan/copycall/%.copycall
	$(call LSCRIPT_MEM,4G,6G,"$(RSCRIPT) $(SEGMENTCNV) --centromereFile=$(CENTROMERE_TABLE2) --prefix=varscan/segment/$* $<")

include ~/share/modules/gatk.mk
