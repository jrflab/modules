# Compare vcf files
##### DEFAULTS ######
include modules/Makefile.inc

LOGDIR = log/vcfComp.$(NOW)
##### MAKE INCLUDES #####

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all variant_eval gt_concordance

FILTER_SUFFIX := dp_ft.target_ft
ifdef TARGETS_FILE
FILTER_SUFFIX := $(FILTER_SUFFIX).target_ft
endif
ifdef NORMAL_VCF
FILTER_SUFFIX := nft.$(FILTER_SUFFIX)
endif
#VARIANT_TYPES := gatk_snps snvmix2
EVAL_TYPES ?= rnaseq_gatk_snps
COMP_TYPES ?= exonseq_museq exonseq_mutect

all : variant_eval

gt_concordance : $(foreach sample,$(SAMPLES),cmp_vcf/grp/$(sample).gt_concord.grp)

variant_eval : $(foreach sample,$(SAMPLES),cmp_vcf/grp/$(sample).variant_eval.grp)

cmp_vcf/grp/%.gt_concord.grp : $(foreach type,$(VARIANT_TYPES),vcf/%.$(type).$(FILTER_SUFFIX).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T GenotypeConcordance -R $(REF_FASTA) $(foreach i,$^,--eval:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i )  $(foreach i,$^,--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@")
	
cmp_vcf/grp/%.variant_eval.grp : $(foreach type,$(EVAL_TYPES),vcf/%.$(type).vcf) $(foreach type,$(COMP_TYPES),vcf/%.$(type).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T VariantEval --dbsnp $(DBSNP) -R $(REF_FASTA) $(foreach i,$(EVAL_TYPES),--eval:$i vcf/$*.$i.vcf ) $(foreach i,$(COMP_TYPES),--comp:$i vcf/$*.$i.vcf ) -o $@")
#$(call LSCRIPT_MEM,4G,6G,"$(call GATK_MEM,4G) -T VariantEval --dbsnp $(DBSNP) -R $(REF_FASTA)  --eval:$(<F:.$(FILTER_SUFFIX).vcf=) $< $(foreach i,$(wordlist 2,$(words $^),$^),--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@")

include modules/variant_callers/gatk.inc
