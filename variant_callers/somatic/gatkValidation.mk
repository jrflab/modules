# use GATK to gather allelic depths at somatic positions
include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/gatk_validation.$(NOW)

SOMATIC_BED ?= somatic.bed

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : val_vcfs

val_vcfs : $(foreach sample,$(SAMPLES),vcf/$(sample).gatkval_snps.vcf)

vcf/%.gatkval_snps.vcf : bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) -I $< -L $(SOMATIC_BED) -o $@ --output_mode EMIT_ALL_SITES")

