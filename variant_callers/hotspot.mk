# run unified genotyper on hotspots

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/hotspot.$(NOW)
PHONY += hotspot hotspot_vcfs hotspot_tables

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

HOTSPOT_GATK_OPTS = --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -R $(REF_FASTA) -stand_call_conf 0


hotspot : hotspot_vcfs hotspot_tables
hotspot_vcfs : $(if $(SAMPLE_PAIRS),\
	$(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).hotspot.vcf),\
	$(foreach sample,$(SAMPLES),vcf_ann/$(sample).hotspot.vcf))
hotspot_tables : $(if $(SAMPLE_PAIRS),alltables/allTN.hotspot.tab.txt)

vcf_ann/%.hotspot.vcf : vcf/%.hotspot.ac_ft.hotspot_int_ann.hotspot_ext_ann.vcf
	$(INIT) cp $< $@

vcf/%.hotspot.vcf : $(foreach i,int ext,hotspot/%.hotspot-$i.vcf.gz hotspot/%.hotspot-$i.vcf.gz.tbi)
	$(call RUN,-c -s 2G -m 3G,"$(BCFTOOLS2) concat -a $(filter %.vcf.gz,$^) > $@.tmp && \
		$(call VERIFY_VCF,$@.tmp,$@)")

define hotspot-vcf-tumor-normal-i
hotspot/$1_$2.hotspot-$3.vcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai 
	$$(call RUN,-c -s 9G -m 12G,"$$(call GATK_MEM2,4G) \
		-T UnifiedGenotyper $$(HOTSPOT_GATK_OPTS) -I $$(<) -I $$(<<) \
		-alleles $$(HOTSPOT_VCF.$3) -L $$(HOTSPOT_VCF.$3) -o $$@.tmp && \
		$$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(if $(SAMPLE_PAIRS),$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach i,int ext,\
		$(eval $(call hotspot-vcf-tumor-normal-i,$(tumor.$(pair)),$(normal.$(pair)),$i)))))

define hotspot-vcf-sample-i
hotspot/$1.hotspot-$2.vcf : bam/$1.bam bam/$1.bam.bai
	$$(call RUN,-c -s 9G -m 12G,"$$(call GATK_MEM2,4G) \
		-T UnifiedGenotyper $$(HOTSPOT_GATK_OPTS) -I $$(<) \
		-alleles $$(HOTSPOT_VCF.$2) -L $$(HOTSPOT_VCF.$2) \
		-o $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach sample,$(SAMPLES),\
	$(foreach i,int ext,\
		$(eval $(call hotspot-vcf-sample-i,$(sample),$i))))

include modules/vcf_tools/vcftools.mk
