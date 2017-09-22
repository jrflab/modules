include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/haplotype_caller.$(NOW)

# covariates for recalibration step
COVARIATES = -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate

# defaults
VARIANT_CALL_THRESHOLD = 30
VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL = 99.0
VARIANT_RECAL_ANNOTATIONS = QD MQRankSum FS MQ ReadPosRankSum
override HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -rf BadCigar -stand_call_conf $(VARIANT_CALL_THRESHOLD) -R $(REF_FASTA) --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000

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

GATK_HARD_FILTER_SNPS ?= true
GATK_POOL_SNP_RECAL ?= false
GATK_SPLIT_CHR ?= true

.PHONY: gatk
.SECONDARY:
.DELETE_ON_ERROR:

gatk : $(foreach type,indels snps,vcf/all.gatk_hc.$(type).vcf)

%.intervals : %.vcf
	$(INIT) sed '/^#/d' $< | awk '{print $$1":"$$2 }' > $@

#### if SPLIT_CHR is set to true, we will split gatk processing by chromosome
ifeq ($(GATK_SPLIT_CHR),true)

define chr-variants
gatk/chr_gvcf/%.$1.variants.g.vcf.gz : bam/%.bam  #bam/%.bai
	$$(call RUN,-c -s 8G -m 12G,"$$(call GATK_MEM2,8G) -T HaplotypeCaller \
	-L $1 -I $$< -o $$@ $$(HAPLOTYPE_CALLER_OPTS)")

gatk/chr_vcf/all.$1.variants.vcf.gz : $$(foreach sample,$$(SAMPLES),gatk/chr_gvcf/$$(sample).$1.variants.g.vcf.gz)
	$$(call RUN,-c -s 20G -m 20G -w 24:00:00,"$$(call GATK_MEM2,16G) -T GenotypeGVCFs -R $$(REF_FASTA) \
		$$(foreach vcf,$$^,--variant $$(vcf)) -o $$@")

endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-variants,$(chr))))

gatk/vcf/all.variants.vcf.gz : $(foreach chr,$(CHROMOSOMES),gatk/chr_vcf/all.$(chr).variants.vcf.gz)
	$(call RUN,-c -s 4G -m 6G,"$(call GATK_MEM2,3G) -T CombineVariants --disable_auto_index_creation_and_locking_when_reading_rods --assumeIdenticalSamples $(foreach i,$^, --variant $i) -R $(REF_FASTA) -o $@")


else #### no splitting by chr ####

define hapcall-vcf
gatk/gvcf/$1.variants.g.vcf.gz : bam/$1.bam
	$$(call RUN,-c -s 8G -m 12G,"$$(call GATK_MEM2,8G) -T HaplotypeCaller \
	-I $$< $$(HAPLOTYPE_CALLER_OPTS) -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hapcall-vcf,$(sample))))

gatk/vcf/all.variants.vcf.gz : $(foreach sample,$(SAMPLES),gatk/gvcf/$(sample).variants.g.vcf.gz)
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM2,8G) -T GenotypeGVCFs --disable_auto_index_creation_and_locking_when_reading_rods \
		-R $(REF_FASTA) $(foreach vcf,$^,--variant $(vcf) ) -o $@")

endif # split by chr

gatk/vcf/all.variants.snps.vcf.gz : gatk/vcf/all.variants.vcf.gz
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM2,8G) -T SelectVariants  -R $(REF_FASTA)  --variant $< -o $@ \
	 -selectType SNP")

gatk/vcf/all.variants.indels.vcf.gz : gatk/vcf/all.variants.vcf.gz
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM2,8G) -T SelectVariants  -R $(REF_FASTA)  --variant $< -o $@ \
	 -selectType INDEL")

ifeq ($(GATK_HARD_FILTER_SNPS),true)
gatk/vcf/all.variants.snps.filtered.vcf.gz : gatk/vcf/all.variants.snps.vcf.gz
	$(call RUN,-c -s 9G -m 12G,"$(call GATK_MEM2,8G) -T VariantFiltration -R $(REF_FASTA) $(SNP_FILTERS) -o $@ \
	--variant $<")
else 
# run variant recal function
gatk/vcf/all.variants.snps.recal.vcf.gz : gatk/vcf/all.variants.snps.vcf.gz
	$(call RUN,-c -n 6 -s 4G -m 5G,"$(call GATK_MEM2,22G) -T VariantRecalibrator \
		-R $(REF_FASTA) -nt 6 \
	-resource:hapmap$(,)known=false$(,)training=true$(,)truth=true$(,)prior=15.0 $(HAPMAP) \
		-resource:omni$(,)known=false$(,)training=true$(,)truth=false$(,)prior=12.0 $(OMNI) \
		-resource:dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=8.0 $(DBSNP) \
		-input $< \
		-recalFile $@ \
		-tranchesFile $(patsubst %.vcf.gz,%.tranches,$@) \
		-rscriptFile $(patsubst %.vcf.gz,%.snps.plots.R,$@)")

# apply variant recal function
gatk/vcf/all.variants.snps.filtered.vcf.gz : gatk/vcf/all.variants.snps.vcf.gz gatk/vcf/all.variants.snps.recal.vcf.gz
	$(call RUN,-c -s 9G -m 12G,"$(call GATK_MEM2,8G) -T ApplyRecalibration  -R $(REF_FASTA) \
		-input $< \
		-recalFile $(<<) \
		--ts_filter_level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL) \
		-tranchesFile $(patsubst %.vcf.gz,%.tranches,$(<<)) \
		-o $@")
endif

# hard filter indels %=sample
gatk/vcf/all.variants.indels.filtered.vcf.gz : gatk/vcf/all.variants.indels.vcf.gz
	$(call RUN,-c -s 9G -m 12G,"$(call GATK_MEM2,8G) -T VariantFiltration -R $(REF_FASTA) $(INDEL_FILTERS) -o $@ \
		--variant $<")

vcf/all.gatk_hc.%.vcf : gatk/vcf/all.variants.%.filtered.vcf.gz
	mkdir -p $(@D); ln $< $@


include modules/bam_tools/processBam.mk
include modules/vcf_tools/vcftools.mk
