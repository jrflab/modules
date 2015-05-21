# use samtools/bcftools to gather allelic depths at somatic positions
include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/samtools_validation.$(NOW)

SOMATIC_BED ?= somatic.bed

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : val_vcfs

val_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).val.vcf)

define val-tumor-normal
vcf/$1_$2.val.vcf : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,7G,8G,"$(SAMTOOLS2) mpileup -l $(SOMATIC_BED) -f $(REF_FASTA) -t DP$$(,)DV$$(,)DPR$$(,)INFO/DPR -g $$^ | $(BCFTOOLS2) call -c - > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call val-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

