# run unified genotyper on hotspots

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/hotspotTN.$(NOW)
PHONY += hotspot hotspot_vcfs hotspot_tables

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: $(PHONY)

HOTSPOT_GATK_OPTS = --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -R $(REF_FASTA) -stand_call_conf 0 -stand_emit_conf 0


hotspot : hotspot_vcfs hotspot_tables
hotspot_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).hotspot.ac_ft.vcf)
hotspot_tables : alltables/allTN.hotspot.ac_ft.tab.txt

vcf/%.hotspot.vcf : $(foreach i,$(HOTSPOT_VCF_SEQ),hotspot/%.hotspot.$i.hotspot_ann.vcf.gz hotspot/%.hotspot.$i.hotspot_ann.vcf.gz.tbi)
	$(call LSCRIPT_MEM,2G,3G,"$(BCFTOOLS2) concat -a $(filter %.vcf.gz,$^) > $@")

define hotspot-vcf-tumor-normal-i
hotspot/$1_$2.hotspot.$3.vcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai 
	$$(call LSCRIPT_MEM,9G,13G,"$$(call GATK_MEM,8G) \
		-T UnifiedGenotyper $$(HOTSPOT_GATK_OPTS) -I $$(<) -I $$(<<) \
		-alleles $$(HOTSPOT_VCF.$3) -L $$(HOTSPOT_VCF.$3) -o $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach i,$(HOTSPOT_VCF_SEQ),\
		$(eval $(call hotspot-vcf-tumor-normal-i,$(tumor.$(pair)),$(normal.$(pair)),$i))))

include modules/vcf_tools/vcftools.mk
