# Run mutect on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_OPTS = --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION)
MUTECT = $(JAVA6) -Xmx11G -jar $(MUTECT_JAR) --analysis_type MuTect $(MUTECT_OPTS)

MUT_FREQ_REPORT = modules/variant_callers/somatic/mutectReport.Rmd

..DUMMY := $(shell mkdir -p version; echo "$(MUTECT) &> version/mutect.txt")

LOGDIR ?= log/mutect.$(NOW)

PHONY += mutect mutect_vcfs mutect_mafs ext_output mut_report
mutect : mutect_vcfs ext_output # mutect_mafs
mutect_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).mutect.vcf)
ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
mut_report : mutect/report/index.html mutect/lowAFreport/index.html mutect/highAFreport/index.html

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)



# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect-tumor-normal-chr
mutect/chr_vcf/$1_$2.$3.mutect%vcf mutect/chr_tables/$1_$2.$3.mutect%txt : bam/$1%bam bam/$2%bam
	$$(MKDIR) mutect/chr_tables mutect/chr_vcf; $$(call LSCRIPT_CHECK_MEM,12G,16G,"$$(MUTECT) --enable_extended_output \
		--intervals $3 --reference_sequence $$(REF_FASTA) --dbsnp $$(DBSNP) --input_file:tumor $$< --input_file:normal\
		$$(word 2,$$^) -vcf mutect/chr_vcf/$1_$2.$3.mutect.vcf --out mutect/chr_tables/$1_$2.$3.mutect.txt")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

# merge variant tables 
define ext-mutect-tumor-normal
mutect/tables/$1.mutect.txt : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_tables/$1.$$(chr).mutect.txt)
	$$(INIT) head -2 $$< > $$@; for table in $$^; do sed '1,2d' $$$$table >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ext-mutect-tumor-normal,$(pair))))

mutect/report/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_NAMED_MEM,mutect_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) $^")

mutect/lowAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_NAMED_MEM,mutect_lowaf_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --lowAF $^")

mutect/highAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
	$(call LSCRIPT_NAMED_MEM,mutect_highaf_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --highAF $^")

# merge variants 
#$$(INIT) grep '^##' $$< > $$@; echo "##PEDIGREE=<Derived=$1,Original=$2>" >> $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
define mutect-tumor-normal
vcf/$1_$2.mutect.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_vcf/$1_$2.$$(chr).mutect.vcf)
	$$(call LSCRIPT_MEM,4G,8G,"grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk

