# Merge strelka and varscan indel results
LOGDIR = log/merge_strelka_varscan_indels.$(NOW)

include modules/Makefile.inc
include modules/config.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : strelka_varscan_merge strelka_varscan_merge_vcfs strelka_varscan_merge_tables

strelka_varscan_merge : strelka_varscan_merge_vcfs strelka_varscan_merge_tables
strelka_varscan_merge_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).strelka_varscan_indels.vcf)
strelka_varscan_merge_tables : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach ext,$(TABLE_EXTENSIONS),tables/$(pair).strelka_varscan_indels.$(ext).txt))

%.vcf.gz : %.vcf
	$(call LSCRIPT,"$(BGZIP) $<")

%.vcf.gz.csi : %.vcf.gz
	$(call LSCRIPT,"$(BCFTOOLS2) index $<")

vcf/%.strelka_varscan_indels.vcf : vcf/%.$(call VCF_SUFFIXES,strelka_indels).vcf.gz vcf/%.$(call VCF_SUFFIXES,strelka_indels).vcf.gz.csi vcf/%.$(call VCF_SUFFIXES,varscan_indels).vcf.gz vcf/%.$(call VCF_SUFFIXES,varscan_indels).vcf.gz.csi
	$(call LSCRIPT_MEM,9G,12G,"dir=`mktemp`; $(BCFTOOLS2) isec -p $$dir $(filter %.vcf.gz,$^) && cp $$dir/0002.vcf $@")

include modules/vcf_tools/vcftools.mk
