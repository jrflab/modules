#### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc


# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect-tumor-normal-chr
mutect/chr_vcf/$1_$2.$3.mutect%vcf mutect/chr_tables/$1_$2.$3.mutect%txt : bam/$1%bam bam/$2%bam
#	$$(MKDIR) mutect/chr_tables mutect/chr_vcf; $$(call LSCRIPT_MEM,12G,16G,"$$(MUTECT) --enable_extended_output --intervals $3 --reference_sequence $$(REF_FASTA) --cosmic $$(COSMIC) --dbsnp $$(DBSNP1PC) --input_file:tumor $$< --input_file:normal $$(word 2,$$^) -vcf mutect/chr_vcf/$1_$2.$3.mutect.vcf --out mutect/chr_tables/$1_$2.$3.mutect.txt")
	$$(MKDIR) mutect/chr_tables mutect/chr_vcf; $$(call LSCRIPT_MEM,12G,16G,"$$(MUTECT) --enable_extended_output --intervals $3 --reference_sequence $$(REF_FASTA) --dbsnp $$(DBSNP1PC) --input_file:tumor $$< --input_file:normal $$(word 2,$$^) -vcf mutect/chr_vcf/$1_$2.$3.mutect.vcf --out mutect/chr_tables/$1_$2.$3.mutect.txt")
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

