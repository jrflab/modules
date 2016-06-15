# run unified genotyper on hotspots

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/hotspot.$(NOW)
PHONY += hotspot hotspot_vcfs hotspot_tables

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

HOTSPOT_GATK_OPTS = --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -R $(REF_FASTA) -stand_call_conf 0 -stand_emit_conf 0


hotspot: hotspot_vcfs hotspot_tables
hotspot_vcfs : $(foreach sample,$(SAMPLES),vcf/$(sample).hotspot.ac_ft.vcf)
hotspot_tables : alltables/all.hotspot.ac_ft.tab.txt

vcf/%.hotspot.vcf : $(foreach i,$(HOTSPOT_VCF_SEQ),hotspot/%.hotspot.$i.hotspot_ann.vcf.gz hotspot/%.hotspot.$i.hotspot_ann.vcf.gz.tbi)
	$(call LSCRIPT_MEM,2G,3G,"$(BCFTOOLS2) concat -a $(filter %.vcf.gz,$^) > $@")

define hotspot-vcf-sample-i
hotspot/$1.hotspot.$2.vcf : bam/$1.bam bam/$1.bam.bai
	$$(call LSCRIPT_MEM,8G,10G,"$$(call GATK_MEM,8G) \
		-T UnifiedGenotyper $$(HOTSPOT_GATK_OPTS) -I $$(<) \
		-alleles $$(HOTSPOT_VCF.$2) -L $$(HOTSPOT_VCF.$2) -o $$@")
endef
$(foreach sample,$(SAMPLES),\
	$(foreach i,$(HOTSPOT_VCF_SEQ),\
		$(eval $(call hotspot-vcf-sample-i,$(sample),$i))))

include modules/vcf_tools/vcftools.mk
