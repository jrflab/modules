# GATK_HARD_FILTER_SNPS = true/false (default: true)
# GATK_POOL_SNP_RECAL = true/false (default: false)
# SPLIT_CHR = true/false (default: true)

ifndef GATK_MK

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/gatk.$(NOW)

# defaults

GATK_HARD_FILTER_SNPS ?= true
GATK_POOL_SNP_RECAL ?= false
SPLIT_CHR ?= true

VARIANT_EVAL_GATK_REPORT = $(RSCRIPT) modules/variant_callers/variantEvalGatkReport.R

VARIANT_CALL_THRESHOLD = 30
VARIANT_EMIT_THRESHOLD = 10
VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL = 99.0
VARIANT_RECAL_ANNOTATIONS = QD MQRankSum FS MQ ReadPosRankSum
HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -rf BadCigar -stand_call_conf $(VARIANT_CALL_THRESHOLD) -R $(REF_FASTA)

INDEL_FILTERS = --filterName 'QD' --filterExpression 'QD < 2.0' --filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < -20.0'  --filterName 'InbreedingCoeff' --filterExpression 'InbreedingCoeff < -0.8'  --filterName 'FS' --filterExpression 'FS > 200.0' --filterName 'DP' --filterExpression 'DP < 5'

MQ_THRESHOLD ?= 40.0
QD_THRESHOLD ?= 2.0
FS_THRESHOLD ?= 60.0
HAP_SCORE_THRESHOLD ?= 13.0
MQ_RANKSUM_THRESHOLD ?= -12.5
READPOS_RANKSUM_THRESHOLD ?= -8.0
SNP_FILTERS := --filterName 'QD' --filterExpression 'QD < $(QD_THRESHOLD)' --filterName 'MQ' --filterExpression 'MQ < $(MQ_THRESHOLD)' --filterName 'FS' --filterExpression 'FS > $(FS_THRESHOLD)' --filterName 'HapScore' --filterExpression 'HaplotypeScore > $(HAP_SCORE_THRESHOLD)' --filterName 'MQRankSum' --filterExpression 'MQRankSum < $(MQ_RANKSUM_THRESHOLD)' --filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < $(READPOS_RANKSUM_THRESHOLD)' --filterName 'Depth' --filterExpression 'DP < 5'
SNP_FILTER_INTRON = false
ifeq ($(SNP_FILTER_INTRON),true)
SNP_FILTERS += --filterName 'SnpEff' --filterExpression \"SNPEFF_EFFECT == 'UPSTREAM' || SNPEFF_EFFECT == 'DOWNSTREAM' || SNPEFF_EFFECT == 'INTERGENIC' || SNPEFF_EFFECT == 'INTRAGENIC' || SNPEFF_EFFECT == 'INTRON'\"
endif

REPORT_STRATIFICATION := Filter

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: gatk_vcfs

gatk_vcfs : $(foreach sample,$(SAMPLES),vcf_ann/$(sample).gatk_snps.vcf vcf_ann/$(sample).gatk_indels.vcf)
gatk_vcfs : $(foreach sample,$(SAMPLES),vcf_ann/$(sample).gatk_snps.vcf vcf_ann/$(sample).gatk_indels.vcf)

vcf_ann/%.gatk_snps.vcf : vcf/%.gatk_snps.pass.eff.vcf
	mkdir -p $(@D); ln -f $< $@
vcf_ann/%.gatk_indels.vcf : vcf/%.gatk_indels.pass.eff.vcf
	mkdir -p $(@D); ln -f $< $@

###### RECIPES #######

%.intervals : %.vcf
	$(INIT) sed '/^#/d' $< | awk '{print $$1":"$$2 }' > $@

#### if SPLIT_CHR is set to true, we will split gatk processing by chromosome
ifeq ($(SPLIT_CHR),true)

## call sample sets
ifdef SAMPLE_SET_PAIRS
define hapcall-vcf-sets-chr
gatk/chr_vcf/$1.$2.variants.vcf : $$(foreach sample,$$(samples.$1),gatk/chr_vcf/$$(sample).$2.variants.intervals) $$(foreach sample,$$(samples.$1),bam/$$(sample).bam bam/$$(sample).bai)
	$$(call RUN,-c -s 9G -m 12G -w 24:00:00,"$$(call GATK_MEM,8G) -T HaplotypeCaller $$(HAPLOTYPE_CALLER_OPTS) \
		$$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) $$(foreach intervals,$$(filter %.intervals,$$^),-L $$(intervals) ) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call hapcall-vcf-sets-chr,$(set),$(chr)))))

define merge-chr-variants-sets
gatk/vcf/$1.variants.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.vcf)
	$$(call RUN,-c -s 4G -m 6G,"$$(call GATK_MEM,3G) -T CombineVariants \
	--assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@")
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call merge-chr-variants-sets,$(set))))
endif # def SAMPLE_SETS

define chr-variants
gatk/chr_vcf/%.$1.variants.vcf : bam/%.bam bam/%.bai
	$$(call RUN,-s 8G -m 12G -w 24:00:00,"$$(call GATK_MEM,8G) -T HaplotypeCaller \
	-L $1 -I $$< -o $$@ \
	$$(HAPLOTYPE_CALLER_OPTS)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-variants,$(chr))))

define merge-chr-variants
gatk/vcf/$1.variants.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.vcf)
	$$(call RUN,-s 4G -m 6G -c,"$$(call GATK_MEM,3G) -T CombineVariants --assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-chr-variants,$(sample))))

else #### no splitting by chr ####

## call sample sets
ifdef SAMPLE_SETS
define hapcall-vcf-sets
gatk/vcf/$1.variants.vcf : $$(foreach sample,$$(samples.$1),gatk/vcf/$$(sample).variants.vcf) $$(foreach sample,$$(samples.$1),bam/$$(sample).bam bam/$$(sample).bai)
	$$(call RUN,-s 9G -m 12G -c,"$$(call GATK_MEM,8G) -T HaplotypeCaller -R $$(REF_FASTA) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) $$(foreach vcf,$$(filter %.vcf,$$^),-L $$(vcf) ) -o $$@")
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call hapcall-vcf-sets,$(set))))
endif

define hapcall-vcf
gatk/vcf/$1.variants.vcf : bam/$1.bam bam/$1.bai
	$$(call RUN,-s 8G -m 12G -c,"$$(call GATK_MEM,8G) -T HaplotypeCaller \
	-I $$< --dbsnp $$(DBSNP) -o $$@  -rf BadCigar \
	-stand_call_conf $$(VARIANT_CALL_THRESHOLD) -stand_emit_conf $$(VARIANT_EMIT_THRESHOLD) -R $$(REF_FASTA)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hapcall-vcf,$(sample))))

endif # split by chr

# select snps % = sample
gatk/vcf/%.variants.snps.vcf : gatk/vcf/%.variants.vcf 
	$(call RUN,-s 8G -m 12G,"$(call GATK_MEM,8G) -T SelectVariants  -R $(REF_FASTA)  --variant $<  -o $@ \
	 -selectType SNP")

# select indels % = indels
gatk/vcf/%.variants.indels.vcf : gatk/vcf/%.variants.vcf 
	$(call RUN,-s 8G -m 12G,"$(call GATK_MEM,8G) -T SelectVariants -R $(REF_FASTA) --variant $<  -o $@ \
	-selectType INDEL")

#$(call VARIANT_RECAL,$@,$^)
define VARIANT_RECAL
$(call RUN,-n 6 -s 4G -m 5G,"$(call GATK_MEM,22G) -T VariantRecalibrator \
	-R $(REF_FASTA) -nt 6 \
-resource:hapmap$(,)known=false$(,)training=true$(,)truth=true$(,)prior=15.0 $(HAPMAP) \
	-resource:omni$(,)known=false$(,)training=true$(,)truth=false$(,)prior=12.0 $(OMNI) \
	-resource:dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=8.0 $(DBSNP) \
	$(foreach i,$(VARIANT_RECAL_ANNOTATIONS), -an $i) \
	$(foreach i,$(filter %.vcf,$2), -input $i) \
	-recalFile $1 -tranchesFile $(basename $1).tranches -rscriptFile $(basename $1).snps.plots.R")
endef


# apply variant recal function
# arguments: vcf, recal file
#$(call APPLY_VARIANT_RECAL,$@,input,recal-file)
define APPLY_VARIANT_RECAL
$(call RUN,-s 9G -m 12G,"$(call GATK_MEM,8G) -T ApplyRecalibration  -R $(REF_FASTA) \
	-input $2 \
	-recalFile $3 \
	--ts_filter_level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL) \
	-tranchesFile $(basename $3).tranches -o $1")
endef

# apply variant recal %=sample
ifeq ($(GATK_HARD_FILTER_SNPS),true)
gatk/vcf/%.variants.snps.filtered.vcf : gatk/vcf/%.variants.snps.vcf
	$(call RUN,-s 9G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) $(SNP_FILTERS) -o $@ \
	--variant $<")
else 

# pool sample vcfs for recalibration
ifeq ($(GATK_POOL_SNP_RECAL),true)
gatk/vcf/samples.snps.recal.vcf : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.vcf)
	$(call VARIANT_RECAL,$@,$^)
define sample-apply-recal
gatk/vcf/$1.variants.snps.filtered.vcf : gatk/vcf/$1.variants.snps.vcf gatk/vcf/samples.snps.recal.vcf 
	$$(call APPLY_VARIANT_RECAL,$$@,$$<,$$(word 2,$$^))
endef
$(foreach sample,$(SAMPLES),$(eval $(call sample-apply-recal,$(sample))))

ifdef SAMPLE_SETS
gatk/vcf/sets.snps.recal.vcf : $(foreach set,$(SAMPLE_SET_PAIRS),gatk/vcf/$(set).variants.snps.vcf )
	$(call VARIANT_RECAL,$@,$^)
define sets-apply-recal
gatk/vcf/$1.variants.snps.filtered.vcf : gatk/vcf/$1.variants.snps.vcf gatk/vcf/sets.snps.recal.vcf
	$$(call APPLY_VARIANT_RECAL,$$@,$$<,$$(word 2,$$^))
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call sets-apply-recal,$(set))))
endif

else 
gatk/vcf/%.variants.snps.recal.vcf : gatk/vcf/%.variants.snps.vcf
	$(call VARIANT_RECAL,$@,$^)
gatk/vcf/%.variants.snps.filtered.vcf : gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.recal.vcf
	$(call APPLY_VARIANT_RECAL,$@,$<,$(word 2,$^))
endif
endif

# hard filter indels %=sample
gatk/vcf/%.variants.indels.filtered.vcf : gatk/vcf/%.variants.indels.vcf
	$(call RUN,-c -s 9G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) $(INDEL_FILTERS) -o $@ \
	--variant $<")

# filter for only novel snps/indels
%.novel.txt : %.txt
	$(INIT) /bin/awk 'NR == 1 || $$4 == "."' $< > $@


vcf/%.gatk_snps.vcf : gatk/vcf/%.variants.snps.filtered.vcf
	$(INIT) ln -f $< $@

vcf/%.gatk_indels.vcf : gatk/vcf/%.variants.indels.filtered.vcf
	$(INIT) ln -f $< $@


reports/%/index.html : reports/%.dp_ft.grp metrics/hs_metrics.txt
	$(call RUN,,"$(VARIANT_EVAL_GATK_REPORT) --metrics $(word 2,$^) --outDir $(@D) $<")


# merge variants 
include modules/bam_tools/processBam.mk
include modules/vcf_tools/vcftools.mk

endif
GATK_MK = true
