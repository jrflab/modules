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

VARIANT_TYPES = mutect_snps mutect_indels
mutect2_vcfs : $(foreach type,$(VARIANT_TYPES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(type).vcf))

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)


# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect2-tumor-normal-chr
mutect2/chr_vcf/$1_$2.$3.mutect_snps_indels.vcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_CHECK_MEM,15G,17G,"$$(MUTECT2) --intervals $3 -I:tumor $$< -I:normal $$(<<) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect2-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

define mutect2-tumor-normal
mutect2/vcf/$1_$2.mutect_snps_indels.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$$(chr).mutect_snps_indels.vcf)
	$$(call LSCRIPT_MEM,4G,8G,"(grep '^#' $$< | sed 's/^##\(.*\),Number=R/##\1,Number=A/'; \
		cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -) > $$@")

vcf/$1_$2.%.vcf : mutect2/vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect2-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

mutect2/vcf/%.mutect_snps.vcf mutect2/vcf/%.mutect_indels.vcf: mutect2/vcf/%.mutect_snps_indels.vcf
	$(call LSCRIPT_MEM,4G,6G,"$(SPLIT_SNPS_INDELS_VCF) -s mutect2/vcf/$*.mutect_snps.vcf -i mutect2/vcf/$*.mutect_indels.vcf $<")

include modules/vcf_tools/vcftools.mk
