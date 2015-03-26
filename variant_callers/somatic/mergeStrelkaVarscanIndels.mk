# Merge strelka and varscan indel results
LOGDIR = log/merge_strelka_varscan_indels.$(NOW)

include modules/variant_callers/somatic/somaticVariantCaller.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : strelka_varscan_merge strelka_varscan_merge_vcfs strelka_varscan_merge_tables

strelka_varscan_merge : strelka_varscan_merge_vcfs strelka_varscan_merge_tables
strelka_varscan_merge_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).strelka_varscan_indels.vcf)
strelka_varscan_merge_tables : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach ext,$(TABLE_EXTENSIONS),tables/$(pair).strelka_varscan_indels.$(ext).txt))

vcf/%.strelka_varscan_indels.vcf : vcf/%.$(call VCF_SUFFIXES,strelka_indels).vcf vcf/%.$(call VCF_SUFFIXES,varscan_indels).vcf
	$(call LSCRIPT_MEM,8G,10G,"$(call GATK_MEM,8G) -T CombineVariants --variant $< --variant $(<<) -o $@ -genotypeMergeOptions UNIQUIFY -R $(REF_FASTA)")

