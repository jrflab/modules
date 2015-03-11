# Run VarScan on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######


REF ?= hg19
LOGDIR = log/varscan.$(NOW)

SPLIT_CHR ?= true

##### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

FP_FILTER = $(PERL) $(HOME)/share/usr/bin/fpfilter.pl
BAM_READCOUNT = $(HOME)/share/usr/bin/bam-readcount

VARSCAN_TO_VCF = $(PERL) scripts/varscanTNtoVcf.pl


MIN_MAP_QUAL ?= 1

MIN_VAR_FREQ ?= 0.05

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: varscan varscan_vcfs varscan_tables

SNP_VCF_EFF_FIELDS += VAF
INDEL_VCF_EFF_FIELDS += VAF

VARIANT_TYPES = varscan_indels varscan_snps
VCFS = $(foreach pair,$(SAMPLE_PAIRS),\
		   $(foreach suff,$(call VCF_SUFFIXES,$(VARIANT_TYPES)), \
			   vcf/$(pair).$(suff).vcf))
TABLES = $(foreach pair,$(SAMPLE_PAIRS),\
			 $(foreach suff,$(call TABLE_SUFFIXES,$(VARIANT_TYPES)),tables/$(pair).$(suff).txt))
TABLES += $(foreach suff,$(call TABLE_SUFFIXES,$(VARIANT_TYPES)),alltables/allTN.$(suff).txt)

varscan : varscan_vcfs varscan_tables
varscan_vcfs : $(VCFS)
varscan_tables : $(TABLES)

%.Somatic.txt : %.txt
	$(call LSCRIPT_MEM,5G,8G,"$(call VARSCAN_MEM,4G) somaticFilter $< && $(call VARSCAN_MEM,4G) processSomatic $< && rename .txt.Somatic .Somatic.txt $** && rename .txt.Germline .Germline.txt $** && rename .txt.LOH .LOH.txt $** && rename .txt.hc .hc.txt $**")

ifeq ($(SPLIT_CHR),true)
define varscan-somatic-tumor-normal-chr
varscan/chr_tables/$1_$2.$3.varscan_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) somatic \
	<($$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$(word 2,$$^)) \
	<($$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--min-var-freq $(MIN_VAR_FREQ) --output-indel varscan/chr_tables/$1_$2.$3.indel.txt --output-snp varscan/chr_tables/$1_$2.$3.snp.txt && touch $$@")

varscan/chr_tables/$1_$2.$3.indel.txt : varscan/chr_tables/$1_$2.$3.varscan_timestamp
varscan/chr_tables/$1_$2.$3.snp.txt : varscan/chr_tables/$1_$2.$3.varscan_timestamp

varscan/chr_tables/$1_$2.$3.%.fp_pass.txt : varscan/chr_tables/$1_$2.$3.%.txt bam/$1.bam
	$$(call LSCRIPT_MEM,8G,35G,"$$(FP_FILTER) --output-basename varscan/chr_tables/$1_$2.$3.$$* $$< <($$(BAM_READCOUNT) -f $$(REF_FASTA) $$(word 2,$$^) $3) && head -1 $$< > $$@ && cat varscan/chr_tables/$1_$2.$3.$$*.pass >> varscan/chr_tables/$1_$2.$3.$$*.fp_pass.txt")

endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach i,$(SETS_SEQ), \
		$(foreach tumor,$(call get_tumors,$(set.$i)), \
			$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor),$(call get_normal,$(set.$i)),$(chr))))))

define merge-varscan-tables
varscan/tables/$1.%.txt : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_tables/$1.$$(chr).%.txt)
	$$(INIT) head -1 $$< > $$@ && for x in $$^; do sed 1d $$$$x >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call merge-varscan-tables,$(pair))))

define convert-varscan-tumor-normal
varscan/vcf/$1_$2.%.vcf : varscan/tables/$1_$2.%.txt
	$$(INIT) $$(VARSCAN_TO_VCF) -f $$(REF_FASTA) -t $1 -n $2 $$< | $$(VCF_SORT) $$(REF_DICT) - > $$@

endef
$(foreach i,$(SETS_SEQ), \
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call convert-varscan-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

else # no splitting by chr

define varscan-somatic-tumor-normal
varscan/tables/$1_$2.varscan_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) somatic \
	<($$(SAMTOOLS) mpileup -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$(word 2,$$^)) \
	<($$(SAMTOOLS) mpileup -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--min-var-freq $(MIN_VAR_FREQ) --output-indel varscan/chr_tables/$1_$2.indel.txt --output-snp varscan/chr_tables/$1_$2.snp.txt && touch $$@")

varscan/tables/$1_$2.indel.txt : varscan/tables/$1_$2.varscan_timestamp
varscan/tables/$1_$2.snp.txt : varscan/tables/$1_$2.varscan_timestamp

varscan/tables/$1_$2.%.fp_pass.txt : varscan/tables/$1_$2.%.txt bamrc/$1.bamrc
	$$(call LSCRIPT_MEM,8G,45G,"$$(FP_FILTER) --output-basename varscan/tables/$1_$2.$$* $$^ && mv varscan/tables/$1_$2.$$*.pass varscan/tables/$1_$2.$$*.fp_pass.txt")

endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call varscan-somatic-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))
endif

vcf/%.varscan_indels.vcf : varscan/vcf/%.indel.Somatic.vcf
	$(INIT) ln $< $@

vcf/%.varscan_snps.vcf : varscan/vcf/%.snp.Somatic.vcf
	$(INIT) ln $< $@

define bamrc-chr
bamrc/%.$1.chr_bamrc : bam/%.bam
	$$(call LSCRIPT_MEM,2G,3G,"$$(BAM_READCOUNT) -f $$(REF_FASTA) $$< $1 > $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call bamrc-chr,$(chr))))

include modules/variant_callers/gatk.mk
