# Compare vcf files
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/vcfComp.$(NOW)
SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all variant_eval gt_concordance


FILTER_SUFFIX := dp_ft
ifdef NORMAL_VCF
FILTER_SUFFIX := nft.$(FILTER_SUFFIX)
endif
VARIANT_TYPES := gatk_snps snvmix2

all : gt_concordance variant_eval

gt_concordance : $(foreach sample,$(SAMPLES),cmp_vcf/grp/$(sample).gt_concord.grp)

variant_eval : $(foreach sample,$(SAMPLES),cmp_vcf/grp/$(sample).variant_eval.grp)

cmp_vcf/grp/%.gt_concord.grp : $(foreach type,$(VARIANT_TYPES),vcf/%.$(type).$(FILTER_SUFFIX).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T GenotypeConcordance -R $(REF_FASTA) $(foreach i,$^,--eval:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i )  $(foreach i,$^,--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@")
	
cmp_vcf/grp/%.variant_eval.grp : $(foreach type,$(VARIANT_TYPES),vcf/%.$(type).$(FILTER_SUFFIX).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T VariantEval --dbsnp $(DBSNP) -R $(REF_FASTA) $(foreach i,$^,--eval:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) $(foreach i,$^,--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@")
#$(call LSCRIPT_MEM,4G,6G,"$(call GATK_MEM,4G) -T VariantEval --dbsnp $(DBSNP) -R $(REF_FASTA)  --eval:$(<F:.$(FILTER_SUFFIX).vcf=) $< $(foreach i,$(wordlist 2,$(words $^),$^),--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@")

include ~/share/modules/gatk.inc
