# vim:set ft=make:
#

ifndef GATK_INC
DEFAULT_JAVA_MEM = 18G
GATK_JAR = $(HOME)/share/usr/lib/java/GenomeAnalysisTK.jar
GATK_JAR2 = $(HOME)/share/usr/lib/java/GenomeAnalysisTK-3.3-0.jar
GATK = $(call GATK_MEM,$(DEFAULT_JAVA_MEM))
GATK_MEM = $(JAVA7) -Xmx$(1) -jar $(GATK_JAR) -S LENIENT
GATK2 = $(call GATK_MEM2,$(DEFAULT_JAVA_MEM))
GATK_MEM2 = $(JAVA7) -Xmx$(1) -jar $(GATK_JAR2) -S LENIENT
INDEL_REALIGN = $(JAVA7) -Xmx$(DEFAULT_JAVA_MEM) -jar $(GATK_JAR) -S LENIENT -T IndelRealigner

#VCF_SAMPLES ?= 0
VCF_GEN_IDS ?= GT AD DP
VCF_FIELDS = CHROM POS ID REF ALT FILTER
# vcf genotype fields
#VCF_FIELDS += $(foreach gen,$(VCF_GEN_IDS),$(foreach i,$(shell seq 0 $$(($1-1))),'GEN[$i].$(gen)'))
# vcf eff fields
ANN_FIELDS = 'ANN[*].ALLELE' 'ANN[*].EFFECT' 'ANN[*].IMPACT' 'ANN[*].GENE' 'ANN[*].GENEID' 'ANN[*].FEATURE' 'ANN[*].FEATUREID' 'ANN[*].BIOTYPE' 'ANN[*].RANK' 'ANN[*].HGVS_C' 'ANN[*].HGVS_P' 'ANN[*].CDNA_POS' 'ANN[*].CDNA_LEN' 'ANN[*].CDS_POS' 'ANN[*].CDS_LEN' 'ANN[*].AA_POS' 'ANN[*].AA_LEN' 'ANN[*].DISTANCE'
#VCF_FIELDS += $(addprefix dbNSFP_,$(NSFP_FIELDS))
#VCF_FIELDS += HRun
# chasm
#VCF_FIELDS += chasm_mut chasm_pval chasm_score chasm_fdr
# fathmm
#VCF_FIELDS += fathmm_query fathmm_pred fathmm_score
# haplotype insufficiency
#VCF_FIELDS += hap_insuf
# transfic
#VCF_FIELDS += transficSIFT_score transficSIFT_pred transficPolyphen2_HVAR_score transficPolyphen2_HVAR_pred transficMutationAssessor_score transficMutationAssessor_pred

#VCF_ANNOTATIONS = RMSMappingQuality QualByDepth MappingQualityZero MappingQualityRankSumTest HaplotypeScore FisherStrand Coverage BaseQualityRankSumTest HomopolymerRun
VCF_ANNOTATIONS = QualByDepth HaplotypeScore FisherStrand Coverage BaseQualityRankSumTest DepthPerAlleleBySample

#$(call ALL_VCF_EFF_FIELDS,NUM_SAMPLES)
SNP_VCF_EFF_FIELDS = $(VCF_FIELDS) $(call VCF_GEN_FIELDS,$1) $(VCF_ANN_FIELDS) $(VCF_EFF_FIELDS) $(VCF_NSFP_FIELDS) $(CHASM_FIELDS) $(FATHMM_FIELDS) $(TRANSFIC_FIELDS)
INDEL_VCF_EFF_FIELDS = $(VCF_FIELDS) $(call VCF_GEN_FIELDS,$1) $(VCF_ANN_FIELDS) $(VCF_EFF_FIELDS) 
#ALL_VCF_GATK_FIELDS = $(VCF_FIELDS) $(call VCF_GEN_FIELDS,$1) $(VCF_GATK_FIELDS) $(VCF_NSFP_FIELDS) $(VCF_ANN_FIELDS)

# my scripts
VCF_TO_TABLE = ~/share/scripts/vcfToTable.R
INTRON_POSN_LOOKUP = ~/share/scripts/posnGeneLookup.pl
RBIND = ~/share/scripts/rbind.R

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

TABLE_FIELDS = CHROM POS ID REF ALT AD DP PL AC AF AN BaseQRankSum DB Dels FS HaplotypeScore MLEAC MLEAF MQ MQ0 MQRankSum QD SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT SNPEFF_GENE_BIOTYPE SNPEFF_GENE_NAME SNPEFF_TRANSCRIPT_ID

GT_TABLE_FIELDS = GT AD DP PL

REPORT_STRATIFICATION := Filter
endif
GATK_INC = true
