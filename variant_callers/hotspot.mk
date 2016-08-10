# run unified genotyper on hotspots

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/hotspot.$(NOW)
PHONY += hotspot hotspot_vcfs hotspot_tables

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

HOTSPOT_GATK_OPTS = --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -R $(REF_FASTA) -stand_call_conf 0 -stand_emit_conf 0


hotspot : hotspot_vcfs # hotspot_tables
hotspot_vcfs : $(if $(SAMPLE_PAIRS),\
	$(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).hotspot.vcf),\
	$(foreach sample,$(SAMPLES),vcf_ann/$(sample).hotspot.vcf))
hotspot_tables : $(if $(SAMPLE_PAIRS),alltables/allTN.hotspot.tab.txt\
	alltables/all.hotspot.tab.txt)

vcf_ann/%.hotspot.vcf : vcf/%.hotspot.ac_ft.hotspot_ann.vcf
	$(INIT) cp $< $@

vcf/%.hotspot.vcf : $(foreach i,$(HOTSPOT_VCF_SEQ),hotspot/%.hotspot.$i.vcf.gz hotspot/%.hotspot.$i.vcf.gz.tbi)
	$(call LSCRIPT_MEM,2G,3G,"$(BCFTOOLS2) concat -a $(filter %.vcf.gz,$^) > $@")

define hotspot-vcf-tumor-normal-i
hotspot/$1_$2.hotspot.$3.vcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai 
	$$(call LSCRIPT_MEM,8G,10G,"$$(call GATK_MEM,8G) \
		-T UnifiedGenotyper $$(HOTSPOT_GATK_OPTS) -I $$(<) -I $$(<<) \
		-alleles $$(HOTSPOT_VCF.$3) -L $$(HOTSPOT_VCF.$3) -o $$@")
endef
$(if $(SAMPLE_PAIRS),$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach i,$(HOTSPOT_VCF_SEQ),\
		$(eval $(call hotspot-vcf-tumor-normal-i,$(tumor.$(pair)),$(normal.$(pair)),$i)))))

define hotspot-vcf-i
hotspot/%.hotspot.$1.vcf : bam/%.bam bam/%.bam.bai
	$$(call LSCRIPT_MEM,8G,10G,"$$(call GATK_MEM,8G) \
		-T UnifiedGenotyper $$(HOTSPOT_GATK_OPTS) -I $$(<) \
		-alleles $$(HOTSPOT_VCF.$1) -L $$(HOTSPOT_VCF.$1) -o $$@")
endef
$(foreach i,$(HOTSPOT_VCF_SEQ),\
	$(eval $(call hotspot-vcf-i,$(sample),$i)))

include modules/vcf_tools/vcftools.mk
