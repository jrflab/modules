# vim:set ft=make:


DEFAULT_JAVA_MEM = 18G
GATK_JAR = /opt/common/gatk/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar
GATK = $(call GATK_MEM,$(DEFAULT_JAVA_MEM))
GATK_MEM = $(JAVA) -Xmx$(1) -jar $(GATK_JAR) -S LENIENT
INDEL_REALIGN = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(GATK_JAR) -S LENIENT -T IndelRealigner


# picard tools
PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
ANALYZE_COV = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/AnalyzeCovariates.jar $(PICARD_OPTS)
SORT_SAM = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/SortSam.jar $(PICARD_OPTS)
SORT_SAM_MEM = $(JAVA) -Xmx$1 -jar $(PICARD_DIR)/SortSam.jar $(PICARD_OPTS)
REORDER_SAM = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/ReorderSam.jar $(PICARD_OPTS)
REORDER_SAM_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/ReorderSam.jar $(PICARD_OPTS)
MARK_DUP = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/MarkDuplicates.jar $(PICARD_OPTS)
MERGE_SAMS = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/MergeSamFiles.jar $(PICARD_OPTS)
INTERVAL_LIST_TOOL = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/IntervalListTools.jar $(PICARD_OPTS)
CREATE_SEQ_DICT = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/CreateSequenceDictionary.jar
ADD_RG = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/AddOrReplaceReadGroups.jar $(PICARD_OPTS)
VALIDATE_SAM = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/ValidateSamFile.jar

# SNP EFF
SNP_EFF = $(call SNP_EFF_MEM,$(DEFAULT_JAVA_MEM))
SNP_EFF_JAR = $(HOME)/share/usr/snpEff/snpEff.jar
SNP_EFF_CONFIG = $(HOME)/share/usr/snpEff/snpEff.config
SNP_EFF_FLAGS ?= -ud 0 -no-intron -no-intergenic
SNP_EFF_MEM = $(JAVA) -Xmx$1 -jar $(SNP_EFF_JAR) eff -v -c $(SNP_EFF_CONFIG) $(SNP_EFF_FLAGS)
SNP_SIFT_JAR = $(HOME)/share/usr/snpEff/SnpSift.jar
SNP_SIFT_MEM = $(JAVA) -Xmx$1 -jar $(SNP_SIFT_JAR)
SNP_SIFT = $(call SNP_SIFT_MEM,$(DEFAULT_JAVA_MEM))
DB_NSFP = $(HOME)/share/reference/snpEff/dbNSFP2.0.txt
VCF_EFF_ONE_PER_LINE = $(HOME)/share/usr/snpEff/scripts/vcfEffOnePerLine.pl
VCF_EXTRACT_FIELDS ?= CHROM POS ID REF ALT FILTER "GEN[0].GT" "GEN[0].AD" "GEN[0].BQ" "GEN[0].DP" "GEN[0].FA" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK" dbnsfpPolyphen2_HVAR_pred dbnsfpInterpro_domain
GATK_VCF_EXTRACT_FIELDS ?= CHROM POS ID REF ALT FILTER "GEN[0].GT" "GEN[0].AD" "GEN[0].BQ" "GEN[0].DP" "GEN[0].FA" dbnsfpPolyphen2_HVAR_pred dbnsfpInterpro_domain SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT SNPEFF_GENE_BIOTYPE SNPEFF_GENE_NAME SNPEFF_TRANSCRIPT_ID

# my scripts
VCF_TO_TABLE = ~/share/scripts/vcfToTable.R
INTRON_POSN_LOOKUP = ~/share/scripts/posnGeneLookup.pl
RBIND = ~/share/scripts/rbind.R

# covariates for recalibration step
COVARIATES = -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate

# defaults
VARIANT_CONF_THRESHOLD = 20
VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL = 99.0
VARIANT_RECAL_ANNOTATIONS = QD HaplotypeScore MQRankSum ReadPosRankSum MQ

INDEL_FILTERS = --filterName "QD" --filterExpression "QD < 2.0" --filterName "ReadPosRankSum" --filterExpression "ReadPosRankSum < -20.0"  --filterName "InbreedingCoeff" --filterExpression "InbreedingCoeff < -0.8"  --filterName "FS" --filterExpression "FS > 200.0" --filterName "DP" --filterExpression "DP < 5"

SNP_FILTERS = --filterName "QD" --filterExpression "QD < 2.0" --filterName "MQ" --filterExpression "MQ < 20.0" --filterName "FS" --filterExpression "FS > 60.0" --filterName "HapScore" --filterExpression "HaplotypeScore > 13.0" --filterName "MQRankSum" --filterExpression "MQRankSum < -12.5" --filterName "ReadPosRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "Depth" --filterExpression "DP < 5"
SNP_FILTER_INTRON = false
ifeq ($(SNP_FILTER_INTRON),true)
SNP_FILTERS += --filterName "SnpEff" --filterExpression 'SNPEFF_EFFECT == "UPSTREAM" || SNPEFF_EFFECT == "DOWNSTREAM" || SNPEFF_EFFECT == "INTERGENIC" || SNPEFF_EFFECT == "INTRAGENIC" || SNPEFF_EFFECT == "INTRON"'
endif

TABLE_FIELDS = CHROM POS ID REF ALT AD DP PL AC AF AN BaseQRankSum DB Dels FS HaplotypeScore MLEAC MLEAF MQ MQ0 MQRankSum QD SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT SNPEFF_GENE_BIOTYPE SNPEFF_GENE_NAME SNPEFF_TRANSCRIPT_ID

GT_TABLE_FIELDS = GT AD DP PL

REPORT_STRATIFICATION := Filter
