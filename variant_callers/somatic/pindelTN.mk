# Run pindel
# Detect indels
##### DEFAULTS ######
LOGDIR ?= log/pindelTN.$(NOW)

##### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

PINDEL = $(HOME)/share/usr/bin/pindel
PINDEL2VCF = $(HOME)/share/usr/bin/pindel2vcf

PINDEL2VCF_OPTS = -G -co 50 -ir 2 -il 3 -pr 2 -pl 3
PINDEL_SOMATIC_AD_FILTER_VCF = python modules/vcf_tools/somatic_ad_filter_vcf.py
PINDEL_SOURCE_ANN_VCF = python modules/vcf_tools/annotate_source_vcf.py --source pindel

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : pindel_vcfs

PINDEL_VCFS = $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).pindel.vcf)
#TABLES += $(foreach suff,$(TABLE_SUFFIXES),tables/all.$(suff).txt)

pindel : pindel_vcfs
pindel_vcfs : $(PINDEL_VCFS)

pindel/ins_size/%.insert_size.txt : bam/%.bam
	$(call RUN,,"$(SAMTOOLS) view $< | head -10000 | $(GET_INSERT_SIZE) - > $@")

define pindel-config-tumor-normal
pindel/config/$1_$2.pindel_config.txt : bam/$1.bam bam/$2.bam pindel/ins_size/$1.insert_size.txt pindel/ins_size/$2.insert_size.txt
	$$(INIT) rm -f $$@; for sample in $1 $2; do echo "bam/$$$$sample.bam " `perl -ne 'if (/^Read span: mean ([0-9]+)/) { print "$$$$1"}' pindel/ins_size/$$$$sample.insert_size.txt` " $$$$sample" >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call pindel-config-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define pindel-chr
pindel/%.$1.pindel_timestamp : pindel/config/%.pindel_config.txt
	$$(call RUN,-n 4 -s 3G -m 6G,"$$(MKDIR) pindel/$$* && $$(PINDEL) -T 4 -f $$(REF_FASTA) -i $$< -o pindel/$$*/$1 -c $1 && touch $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call pindel-chr,$(chr))))

pindel/%.pindel_timestamp : $(foreach chr,$(CHROMOSOMES),pindel/%.$(chr).pindel_timestamp)
	$(INIT) touch $@

define pindel-chr-vcf
pindel/chr_vcf/%.pindel_$1.vcf : pindel/%.$1.pindel_timestamp
	$$(call RUN,-c -s 3G -m 4G,"$$(PINDEL2VCF) -P pindel/$$*/$1 -v $$@.tmp -c $1 -r $$(REF_FASTA) -R $$(REF_NAME) -d $$(REF_DATE) $$(PINDEL2VCF_OPTS) && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call pindel-chr-vcf,$(chr))))

define merge-pindel-chr
pindel/vcf/$1.pindel.vcf : $$(foreach chr,$$(CHROMOSOMES),pindel/chr_vcf/$1.pindel_$$(chr).vcf)
	$$(call RUN,-s 4G -m 6G,"(grep '^#' $$<; grep -v '^#' $$^ ) | $$(VCF_SORT) $$(REF_DICT) - > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call merge-pindel-chr,$(pair))))

define pindel-somatic-filter
vcf/$1_$2.pindel.vcf : pindel/vcf/$1_$2.pindel.vcf
	$$(call RUN,-c -s 4G -m 7G,"$$(PINDEL_SOMATIC_AD_FILTER_VCF) --tumor $1 --normal $2 --pass_only $$< | $$(PINDEL_SOURCE_ANN_VCF) > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call pindel-somatic-filter,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk
