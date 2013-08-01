# Compare vcf files
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/vcfCompTN.$(NOW)
SAMPLE_PAIR_FILE ?= sample_pairs.txt
SAMPLE_FILE ?= samples.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

#SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))


##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all variant_eval gt_concordance


FILTER_SUFFIX := dp_ft
VARIANT_TYPES := mutect museq


all : gt_concordance variant_eval

gt_concordance : $(foreach sample,$(TUMOR_SAMPLES),cmp_vcf/grp/$(sample).gt_concord.grp)

variant_eval : $(foreach sample,$(TUMOR_SAMPLES),cmp_vcf/grp/$(sample).variant_eval.grp)


define select-tumor-variants
cmp_vcf/vcf/$1.%.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_MEM,4G,6G,"$$(call GATK_MEM,3G) -T SelectVariants \
	-R $$(REF_FASTA)  --variant $$<  -o $$@ -sn $1 &> $$(LOG)")
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call select-tumor-variants,$(tumor),$(normal_lookup.$(tumor)))))

cmp_vcf/grp/%.gt_concord.grp : $(foreach type,$(VARIANT_TYPES),cmp_vcf/vcf/%.$(type).$(FILTER_SUFFIX).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T GenotypeConcordance --dbsnp $(DBSNP) -R $(REF_FASTA)  --eval:$(<F:.$(FILTER_SUFFIX).vcf=) $< $(foreach i,$(wordlist 2,$(words $^),$^),--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@ &> $(LOG)")
	
cmp_vcf/grp/%.variant_eval.grp : $(foreach type,$(VARIANT_TYPES),cmp_vcf/vcf/%.$(type).$(FILTER_SUFFIX).vcf)
	$(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T VariantEval --dbsnp $(DBSNP) -R $(REF_FASTA) $(foreach i,$^,--eval:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) $(foreach i,$^,--comp:$(notdir $(i:.$(FILTER_SUFFIX).vcf=)) $i ) -o $@ &> $(LOG)")

include ~/share/modules/gatk.inc
