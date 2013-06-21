# common GATK steps for variant calling
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
#
# GATK needs sorted, properly ordered bams
# INPUTS: sample bams
# OUTPUTS: vcf file for SNPs and indels
# OPTIONS: HARD_FILTER_SNPS = true/false (default: true)
# 		   POOL_SNP_RECAL = true/false (default: false)
# 		   SPLIT_CHR = true/false (default: true)
#

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

HARD_FILTER_SNPS ?= true
POOL_SNP_RECAL ?= false
SPLIT_CHR ?= true

#SNP_EFF_FLAGS ?= -ud 0 -no-intron -no-intergenic

###### RECIPES #######

#### if SPLIT_CHR is set to true, we will split gatk processing by chromosome
ifeq ($(SPLIT_CHR),true)
# $(eval $(call chr-variants,chromosome))
define chr-variants
gatk/chr_vcf/%.$1.variants.vcf gatk/chr_vcf/%.$1.variants.vcf.idx : bam/%.bam bam/%.bai
	$$(call INIT_PARALLEL_MEM,2,5G,6G) $$(call GATK_MEM,8G) -T UnifiedGenotyper \
	-nt 2 -I $$< --dbsnp $$(DBSNP) -o gatk/chr_vcf/$$*.$1.variants.vcf  -glm BOTH -rf BadCigar \
	-stand_call_conf $$(VARIANT_CONF_THRESHOLD)  -R $$(REF_FASTA) &> $$(LOGDIR)/$$@.log
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-variants,$(chr))))

# merge variants 
define merge-chr-variants
gatk/vcf/$1.variants%vcf gatk/vcf/$1.variants%vcf.idx : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.vcf)
	$$(call INIT_MEM,4G,6G) $$(call GATK_MEM,3G) -T CombineVariants --assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o gatk/vcf/$1.variants.vcf &> $$(LOG)
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-chr-variants,$(sample))))

#$(call INIT_MEM,2G,3G) grep '^#' $< > $@; cat $^ | grep -v '^#' | vcfsorter.pl $(REF_DICT) - >> $@ 2> $(LOGDIR)/$@.log

else #### no splitting by chr ####
gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx : bam/%.bam bam/%.bai
	$(call INIT_PARALLEL_MEM,2,4G,6G) $(call GATK_MEM,8G) -T UnifiedGenotyper \
	-nt 2 -I $< --dbsnp $(DBSNP) -o gatk/vcf/$*.variants.vcf  -glm BOTH -rf BadCigar \
	-stand_call_conf $(VARIANT_CONF_THRESHOLD)  -R $(REF_FASTA) &> $(LOG)
endif

# select snps % = sample
gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.vcf.idx : gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	$(call INIT_MEM,8G,12G) $(call GATK_MEM,8G) -T SelectVariants  -R $(REF_FASTA)  --variant $<  -o gatk/vcf/$*.snps.vcf \
	 -selectType SNP &> $(LOGDIR)/$@.log

# select indels % = indels
gatk/vcf/%.variants.indels.vcf gatk/vcf/%.variants.indels.vcf.idx: gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	$(call INIT_MEM,8G,12G) $(call GATK_MEM,8G) -T SelectVariants -R $(REF_FASTA) --variant $<  -o gatk/vcf/$*.indels.vcf \
	-selectType INDEL &> $(LOG)

%.bai : %.bam
	$(call INIT_MEM,4G,8G) $(SAMTOOLS) index $< $@

$(REF_FASTA).fai : $(REF_FASTA)
	$(call INIT_MEM,4G,8G) $(SAMTOOLS) faidx $< || ($(RM) $@ && false)

$(REF_FASTA:.fasta=.dict) : $(REF_FASTA)
	$(call INIT_MEM,5G,6G) $(CREATE_SEQ_DICT) REFERENCE=$< OUTPUT=$@

# variant recalibrator for snps
# not recommended for indels (not enough data)
#ifeq ($(POOL_SNP_RECAL),true)
#gatk/variant_recal/%snps.recal gatk/variant_recal/%snps.tranches gatk/variant_recal/%snps.plots.R : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).%snps.annotated.vcf)
#SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,2G,4G) -pe $(PARALLEL_ENV) 6" $(MKDIR) $(@D); \
#$(GATK) -T VariantRecalibrator -R $(REF_FASTA) $(foreach i,$^, -input $i)  -nt 6 \
#-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $(HAPMAP) \
#-resource:omni,known=false,training=true,truth=false,prior=12.0 $(OMNI) \
#-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $(DBSNP)  $(foreach i,$(VARIANT_RECAL_ANNOTATIONS), -an $i) \
#-recalFile gatk/variant_recal/$*snps.recal \
#-tranchesFile gatk/variant_recal/$*snps.tranches \
#-rscriptFile gatk/variant_recal/$*snps.plots.R &> $(LOGDIR)/$@.log || ($(RM) $@ && false)
#else 
#endif

#$(call VARIANT_RECAL,$@,$^)
define VARIANT_RECAL
$$(call INIT_PARALLEL_MEM,6,2.5G,3G) $$(call GATK_MEM,8G) -T VariantRecalibrator \
	-R $$(REF_FASTA) -nt 6 \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $$(HAPMAP) \
	-resource:omni,known=false,training=true,truth=false,prior=12.0 $$(OMNI) \
	-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $$(DBSNP) \
	$$(foreach i,$$(VARIANT_RECAL_ANNOTATIONS), -an $$i) \
	$$(foreach i,$$(filter %.vcf,$2), -input $i) \
	-recalFile $1 -tranchesFile $$(basename $1).tranches -rscriptFile $$(basename $1).snps.plots.R &> $$(LOG)
endef

# apply variant recal function
# arguments: vcf, recal file
#$(call APPLY_VARIANT_RECAL,$@,input,recal-file)
define APPLY_VARIANT_RECAL
$$(call INIT_MEM,9G,12G) $$(call GATK_MEM,8G) -T ApplyRecalibration  -R $$(REF_FASTA) \
	-input $2 \
	-recalFile $3 \
	--ts_filter_level $$(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL) \
	-tranchesFile $$(basename $3).tranches -o $1 &> $$(LOG)
endef

# apply variant recal %=sample
ifeq ($(HARD_FILTER_SNPS),true)
gatk/vcf/%.snps.filtered.vcf gatk/vcf/%.snps.filtered.vcf.idx : gatk/vcf/%.snps.vcf gatk/vcf/%.snps.vcf.idx
	$(call INIT_MEM,9G,12G) \
	$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) $(SNP_FILTERS) -o gatk/vcf/$*.snps.filtered.vcf \
	--variant $< &> $(LOG)
else 
ifeq ($(POOL_SNP_RECAL),true)
gatk/vcf/snps.recal%vcf gatk/vcf/snps.recal.vcf%idx : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.annotated%vcf) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.annotated%vcf.idx)
	$(call VARIANT_RECAL,gatk/vcf/snps.recal.vcf,$^)
gatk/vcf/%.variants.snps.filtered.vcf gatk/vcf/%.variants.snps.filtered.vcf.idx: gatk/vcf/%.variants.snps.vcf gatk/vcf/snps.recal.vcf gatk/vcf/snps.recal.vcf.idx gatk/vcf/%.snps.vcf.idx 
	$(call APPLY_VARIANT_RECAL,gatk/vcf/$*.variants.snps.filtered.vcf,$<,$(word 2,$^))
#gatk/variant_recal/no_realn.snps.recal : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).snps.annotated.vcf)
#$(VARIANT_RECAL)
#gatk/variant_recal/%.snps.annotated.filtered.vcf : gatk/vcf/%.snps.annotated.vcf gatk/variant_recal/no_realn.snps.recal
#$(APPLY_VARIANT_RECAL)
else 
gatk/vcf/%.snps.recal.vcf gatk/vcf/%.snps.recal.vcf.idx : gatk/vcf/%.snps.vcf gatk/vcf/%.snps.vcf.idx
	$(call VARIANT_RECAL,gatk/vcf/$*.snps.recal.vcf,$^)
gatk/vcf/%.snps.filtered.vcf gatk/vcf/%.snps.filtered.vcf.idx : gatk/vcf/%.snps.vcf gatk/vcf/%.snps.recal.vcf gatk/vcf/%.snps.vcf.idx gatk/vcf/%.snps.recal.vcf.idx
	$(call APPLY_VARIANT_RECAL,gatk/vcf/$*.snps.filtered.vcf,$<,$(word 2,$^))
endif
endif

# hard filter indels %=sample
gatk/vcf/%.variants.indels.filtered.vcf gatk/vcf/%.variants.indels.filtered.vcf.idx : gatk/vcf/%.indels.vcf gatk/vcf/%.indels.vcf.idx
	$(call INIT_MEM,9G,12G) $(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) $(INDEL_FILTERS) -o gatk/vcf/$*.variants.indels.filtered.vcf \
	--variant $< &> $(LOGDIR)/$@.log || ($(RM) $@ && false)

# filter for only novel snps/indels
%.novel.txt : %.txt
	$(INIT) /bin/awk 'NR == 1 || $$4 == "."' $< > $@

vcf/%.gatk_snps.vcf : gatk/vcf/%.variants.snps.filtered.vcf
	$(INIT) cp $< $@

vcf/%.gatk_indels.vcf : gatk/vcf/%.variants.indels.filtered.vcf
	$(INIT) cp $< $@


include ~/share/modules/processBam.mk
include ~/share/modules/vcftools.mk
