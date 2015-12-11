# Run VarScan on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

LOGDIR = log/varscan.$(NOW)

##### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc
include modules/mut_sigs/mutSigReport.mk

IGNORE_FP_FILTER ?= true

FP_FILTER = $(PERL) $(HOME)/share/usr/bin/fpfilter.pl
BAM_READCOUNT = $(HOME)/share/usr/bin/bam-readcount

VARSCAN_TO_VCF = $(PERL) modules/variant_callers/somatic/varscanTNtoVcf.pl

MIN_MAP_QUAL ?= 1
VALIDATION ?= false
MIN_VAR_FREQ ?= $(if $(findstring false,$(VALIDATION)),0.05,0.000001)

VARSCAN_OPTS = $(if $(findstring true,$(VALIDATION)),--validation 1 --strand-filter 0) --min-var-freq $(MIN_VAR_FREQ)

VPATH ?= bam

VARIANT_TYPES = varscan_indels varscan_snps

PHONY += varscan varscan_vcfs varscan_tables
varscan : varscan_vcfs varscan_tables $(if $(findstring false,$(VALIDATION)),varscan_snps_mutsig_report)
varscan_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)))
varscan_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

$(eval $(call mutsig-report-name-vcfs,varscan_snps,$(call SOMATIC_VCFS,varscan_snps)))


%.Somatic.txt : %.txt
	$(call LSCRIPT_MEM,5G,8G,"$(call VARSCAN_MEM,4G) somaticFilter $< && $(call VARSCAN_MEM,4G) processSomatic $< && rename .txt.Somatic .Somatic.txt $** && rename .txt.Germline .Germline.txt $** && rename .txt.LOH .LOH.txt $** && rename .txt.hc .hc.txt $**")

define varscan-somatic-tumor-normal-chr
varscan/chr_tables/$1_$2.$3.varscan_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) somatic \
	<($$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$(word 2,$$^)) \
	<($$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	$$(VARSCAN_OPTS) \
	--output-indel varscan/chr_tables/$1_$2.$3.indel.txt --output-snp varscan/chr_tables/$1_$2.$3.snp.txt && touch $$@")
varscan/chr_tables/$1_$2.$3.indel.txt : varscan/chr_tables/$1_$2.$3.varscan_timestamp
varscan/chr_tables/$1_$2.$3.snp.txt : varscan/chr_tables/$1_$2.$3.varscan_timestamp

varscan/chr_tables/$1_$2.$3.%.fp_pass.txt : varscan/chr_tables/$1_$2.$3.%.txt bamrc/$1.$3.bamrc.gz
	$$(call LSCRIPT_MEM,8G,55G,"$$(VARSCAN) fpfilter $$< <(zcat $$(<<)) --output-file $$@")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

define merge-varscan-pair-type
varscan/tables/$1.$2.txt : $$(foreach chr,$$(CHROMOSOMES),\
	$$(if $$(findstring true,$$(VALIDATION) $$(IGNORE_FP_FILTER)),\
	varscan/chr_tables/$1.$$(chr).$2.txt,\
	varscan/chr_tables/$1.$$(chr).$2.fp_pass.txt))
	$$(INIT) head -1 $$< > $$@ && for x in $$^; do sed 1d $$$$x >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach type,snp indel,$(eval $(call merge-varscan-pair-type,$(pair),$(type)))))

define convert-varscan-tumor-normal
varscan/vcf/$1_$2.%.vcf : varscan/tables/$1_$2.%.txt
	$$(INIT) $$(VARSCAN_TO_VCF) -f $$(REF_FASTA) -t $1 -n $2 $$< | $$(VCF_SORT) $$(REF_DICT) - > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call convert-varscan-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

vcf/%.varscan_indels.vcf : varscan/vcf/%.indel.Somatic.vcf
	$(INIT) ln -f $< $@

vcf/%.varscan_snps.vcf : varscan/vcf/%.snp.Somatic.vcf
	$(INIT) ln -f $< $@

define bamrc-chr
bamrc/%.$1.bamrc.gz : bam/%.bam
	$$(call LSCRIPT_MEM,8G,12G,"$$(BAM_READCOUNT) -f $$(REF_FASTA) $$< $1 | gzip > $$@ 2> /dev/null")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call bamrc-chr,$(chr))))

include modules/variant_callers/gatk.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: $(PHONY)

