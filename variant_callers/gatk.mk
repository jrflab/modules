
# common GATK steps for variant calling
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
#
# GATK needs sorted, properly ordered bams
# INPUTS: sample bams
# OUTPUTS: vcf file for SNPs and indels
# OPTIONS: GATK_HARD_FILTER_SNPS = true/false (default: true)
# 		   GATK_POOL_SNP_RECAL = true/false (default: false)
# 		   SPLIT_CHR = true/false (default: true)
#

ifndef GATK_MK

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

GATK_HARD_FILTER_SNPS ?= true
GATK_POOL_SNP_RECAL ?= false
GATK_SPLIT_CHR ?= true

###### RECIPES #######

%.intervals : %.vcf
	$(INIT) sed '/^#/d' $< | awk '{print $$1":"$$2 }' > $@

#### if SPLIT_CHR is set to true, we will split gatk processing by chromosome
ifeq ($(GATK_SPLIT_CHR),true)

define chr-variants
gatk/chr_gvcf/%.$1.variants.g.vcf : bam/%.bam  #bam/%.bai
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM2,8G) -T HaplotypeCaller \
	-L $1 -I $$< -o $$@.tmp \
	$$(HAPLOTYPE_CALLER_OPTS) && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-variants,$(chr))))

define merge-chr-variants
gatk/gvcf/$1.variants.g.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.g.vcf)
	$$(call LSCRIPT_CHECK_MEM,4G,6G,"$$(call GATK_MEM2,3G) -T CombineVariants --assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-chr-variants,$(sample))))

else #### no splitting by chr ####

define hapcall-vcf
gatk/gvcf/$1.variants.g.vcf : bam/$1.bam
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM2,8G) -T HaplotypeCaller \
	-I $$< $$(HAPLOTYPE_CALLER_OPTS) -o $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hapcall-vcf,$(sample))))

endif # split by chr

gatk/gvcf/all.variants.g.vcf : $(foreach sample,$(SAMPLES),gatk/gvcf/$(sample).variants.g.vcf)
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM2,8G) -T CombineGVCFs -R $(REF_FASTA) \
		$(foreach vcf,$^,--variant $(vcf)) -o $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

gatk/gvcf/all.variants.snps.g.vcf : gatk/gvcf/all.variants.g.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM2,8G) -T SelectVariants  -R $(REF_FASTA)  --variant $< -o $@.tmp \
	 -selectType SNP && $(call VERIFY_VCF,$@.tmp,$@)")

gatk/gvcf/all.variants.indels.g.vcf : gatk/gvcf/all.variants.g.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM2,8G) -T SelectVariants  -R $(REF_FASTA)  --variant $< -o $@.tmp \
	 -selectType INDEL && $(call VERIFY_VCF,$@.tmp,$@)")

ifeq ($(GATK_HARD_FILTER_SNPS),true)
gatk/gvcf/all.variants.snps.filtered.g.vcf : gatk/gvcf/all.variants.snps.g.vcf
	$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call GATK_MEM2,8G) -T VariantFiltration -R $(REF_FASTA) $(SNP_FILTERS) -o $@ \
	--variant $<")
else 
# run variant recal function
gatk/gvcf/all.variants.snps.recal.g.vcf : gatk/gvcf/all.variants.snps.g.vcf
	$(call LSCRIPT_CHECK_PARALLEL_MEM,6,4G,5G,"$(call GATK_MEM2,22G) -T VariantRecalibrator \
		-R $(REF_FASTA) -nt 6 \
	-resource:hapmap$(,)known=false$(,)training=true$(,)truth=true$(,)prior=15.0 $(HAPMAP) \
		-resource:omni$(,)known=false$(,)training=true$(,)truth=false$(,)prior=12.0 $(OMNI) \
		-resource:dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=8.0 $(DBSNP) \
		-input $< \
		-recalFile $@ \
		-tranchesFile $(patsubst %.g.vcf,%.tranches $@) \
		-rscriptFile $(patsubst %.g.vcf,%.snps.plots.R $@)")

# apply variant recal function
gatk/vcf/all.variants.snps.filtered.g.vcf : gatk/vcf/all.variants.snps.g.vcf gatk/gvcf/all.variants.snps.recal.g.vcf
$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call GATK_MEM2,8G) -T ApplyRecalibration  -R $(REF_FASTA) \
	-input $< \
	-recalFile $(<<) \
	--ts_filter_level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL) \
	-tranchesFile $(patsubst %.g.vcf,%.tranches $(<<)) \
	-o $@")
endif

# hard filter indels %=sample
gatk/gvcf/all.variants.indels.filtered.vcf : gatk/vcf/all.variants.indels.vcf
	$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call GATK_MEM2,8G) -T VariantFiltration -R $(REF_FASTA) $(INDEL_FILTERS) -o $@ \
	--variant $<")

# filter for only novel snps/indels
%.novel.txt : %.txt
	$(INIT) /bin/awk 'NR == 1 || $$4 == "."' $< > $@


vcf/all.gatk_hc.g.vcf : gatk/gvcf/all.variants.snps.filtered.g.vcf gatk/gvcf/all.variants.indels.filtered.g.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM2,8G) -T CombineGVCFs -R $(REF_FASTA) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

# merge variants 
include modules/bam_tools/processBam.mk
include modules/vcf_tools/vcftools.mk

endif
GATK_MK = true
