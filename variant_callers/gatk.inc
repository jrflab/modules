# vim:set ft=make:
#

ifndef GATK_INC
DEFAULT_JAVA_MEM = 18G
GATK = $(call GATK_MEM,$(DEFAULT_JAVA_MEM))
GATK_MEM = $(JAVA7) -Xmx$(1) -jar $(GATK_JAR) -S LENIENT
GATK2 = $(call GATK_MEM2,$(DEFAULT_JAVA_MEM))
GATK_MEM2 = $(JAVA7) -Xmx$(1) -jar $(GATK_JAR2) -S LENIENT
INDEL_REALIGN = $(JAVA7) -Xmx$(DEFAULT_JAVA_MEM) -jar $(GATK_JAR) -S LENIENT -T IndelRealigner

# covariates for recalibration step
COVARIATES = -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate

# defaults
VARIANT_CALL_THRESHOLD = 50
VARIANT_EMIT_THRESHOLD = 10
VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL = 99.0
VARIANT_RECAL_ANNOTATIONS = QD MQRankSum FS MQ ReadPosRankSum
override HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -rf BadCigar -stand_call_conf $(VARIANT_CALL_THRESHOLD) -stand_emit_conf $(VARIANT_EMIT_THRESHOLD) -R $(REF_FASTA)

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
endif
GATK_INC = true
