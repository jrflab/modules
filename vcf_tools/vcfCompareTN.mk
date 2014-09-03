# Compare vcf files
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/vcfCompTN.$(NOW)

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all variant_eval gt_concordance


FILTER_SUFFIX := dp_ft
VARIANT_TYPES := mutect museq

VCF_TYPES = varscan_snps strelka_snps
VCF_SUFFIX.gatk_snps := gatk_snps.dp_ft.som_ft.pass.dbsnp.nsfp.eff.chasm.fathmm.transfic
VCF_SUFFIX.gatk_indels := gatk_indels.dp_ft.som_ft.pass.dbsnp.eff
VCF_SUFFIX.strelka_snps := strelka_snps.pass.dbsnp.nsfp.eff.chasm.fathmm.transfic
VCF_SUFFIX.strelka_indels := strelka_indels.pass.dbsnp.eff
VCF_SUFFIX.varscan_snps := varscan_snps.dp_ft.som_ad_ft.pass.dbsnp.eff.nsfp.chasm.fathmm.transfic
VCF_SUFFIX.varscan_indels := varscan_indels.dp_ft.som_ad_ft.pass.dbsnp.eff
VCF_SUFFIX.mutect := mutect.som_ad_ft.pass.dbsnp.nsfp.eff.chasm.fathmm.transfic
VCF_SUFFIX.som_sniper := som_sniper.ss_dp_ft.ss_ft.rn.som_ad_ft.pass.dbsnp.nsfp.eff.chasm.fathmm.transfic
VCF_SUFFIX.scalpel := scalpel.dbsnp.eff
VCF_SUFFIXES := $(foreach type,$(VCF_TYPES),$(VCF_SUFFIX.$(type)))


all : gt_concordance variant_eval

gt_concordance : $(foreach sample,$(TUMOR_SAMPLES),cmp_vcf/grp/$(sample).gt_concord.grp)

variant_eval : $(foreach sample,$(TUMOR_SAMPLES),cmp_vcf/grp/$(sample).variant_eval.grp)


define select-tumor-variants
cmp_vcf/vcf/$1.%.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_MEM,4G,6G,"$$(call GATK_MEM,3G) -T SelectVariants \
	-R $$(REF_FASTA)  --variant $$<  -o $$@ -sn $1 &> $$(LOG)")
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call select-tumor-variants,$(tumor),$(normal_lookup.$(tumor)))))

cmp_vcf/grp/%.gt_concord.grp : $(foreach type,$(VCF_TYPES),cmp_vcf/vcf/%.$(VCF_SUFFIX.$(type)).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T GenotypeConcordance --dbsnp $(DBSNP) -R $(REF_FASTA)  --eval:$(<F:.$(FILTER_SUFFIX).vcf=) $< $(foreach i,$(wordlist 2,$(words $^),$^),--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@")
	
cmp_vcf/grp/%.variant_eval.grp : $(foreach type,$(VCF_TYPES),cmp_vcf/vcf/%.$(VCF_SUFFIX.$(type)).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T VariantEval --dbsnp $(DBSNP) -R $(REF_FASTA) $(foreach i,$^,--eval:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) $(foreach i,$^,--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@")

include ~/share/modules/variant_callers/gatk.inc
