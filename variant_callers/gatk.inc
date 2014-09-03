# vim:set ft=make:

GATK_INCLUDE = true

DEFAULT_JAVA_MEM = 18G
GATK_JAR = $(HOME)/share/usr/lib/java/GenomeAnalysisTK.jar
GATK = $(call GATK_MEM,$(DEFAULT_JAVA_MEM))
GATK_MEM = $(JAVA7) -Xmx$(1) -jar $(GATK_JAR) -S LENIENT
INDEL_REALIGN = $(JAVA7) -Xmx$(DEFAULT_JAVA_MEM) -jar $(GATK_JAR) -S LENIENT -T IndelRealigner


# picard tools

PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT # MAX_RECORDS_IN_RAM=4000000
ANALYZE_COV = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/AnalyzeCovariates.jar $(PICARD_OPTS)
SORT_SAM = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(PICARD_DIR)/SortSam.jar $(PICARD_OPTS)

SORT_SAM_MEM = $(JAVA) -Xmx$1 -jar $(PICARD_DIR)/SortSam.jar $(PICARD_OPTS) MAX_RECORDS_IN_RAM=$2
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
SNP_EFF_MEM = $(JAVA7) -Xmx$1 -jar $(SNP_EFF_JAR) eff -c $(SNP_EFF_CONFIG)
SNP_SIFT_JAR = $(HOME)/share/usr/snpEff/SnpSift.jar
SNP_SIFT_MEM = $(JAVA7) -Xmx$1 -jar $(SNP_SIFT_JAR)
SNP_SIFT = $(call SNP_SIFT_MEM,$(DEFAULT_JAVA_MEM))
DB_NSFP = $(HOME)/share/reference/snpEff/dbNSFP2.4.txt.gz
VCF_EFF_ONE_PER_LINE = $(HOME)/share/usr/snpEff/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = $(HOME)/share/scripts/joinEff.pl


COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)

#VCF_SAMPLES ?= 0
VCF_GEN_IDS ?= GT AD DP
VCF_FIELDS := CHROM POS ID REF ALT FILTER GMAF
VCF_GEN_FIELDS = $(foreach gen,$(VCF_GEN_IDS),$(foreach i,$(shell seq 0 $$(($1-1))),'GEN[$i].$(gen)'))
VCF_EFF_FIELDS = 'EFF[*].EFFECT' 'EFF[*].IMPACT' 'EFF[*].FUNCLASS' 'EFF[*].CODON' 'EFF[*].AA' 'EFF[*].GENE' 'EFF[*].BIOTYPE' 'EFF[*].CODING' 'EFF[*].TRID' 'EFF[*].RANK'
#VCF_ANNOTATIONS = RMSMappingQuality QualByDepth MappingQualityZero MappingQualityRankSumTest HaplotypeScore FisherStrand Coverage BaseQualityRankSumTest HomopolymerRun
VCF_ANNOTATIONS = QualByDepth HaplotypeScore FisherStrand Coverage BaseQualityRankSumTest DepthPerAlleleBySample
NSFP_FIELDS = Ensembl_geneid Ensembl_transcriptid Uniprot_acc Uniprot_id Interpro_domain SIFT_score Polyphen2_HVAR_score Polyphen2_HVAR_pred 1000Gp1_AF ESP6500_AA_AF ESP6500_EA_AF MutationAssessor_pred MutationAssessor_score FATHMM_score MutationTaster_pred MutationTaster_score
VCF_NSFP_FIELDS = $(addprefix dbNSFP_,$(NSFP_FIELDS))
VCF_ANN_FIELDS = HRun
CHASM_FIELDS = chasm_mut chasm_pval chasm_score chasm_fdr
FATHMM_FIELDS = fathmm_query fathmm_pred fathmm_score
TRANSFIC_FIELDS = transficSIFT_score transficSIFT_pred transficPolyphen2_HVAR_score transficPolyphen2_HVAR_pred transficMutationAssessor_score transficMutationAssessor_pred
VCF_GATK_FIELDS = SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT SNPEFF_GENE_BIOTYPE SNPEFF_GENE_NAME SNPEFF_TRANSCRIPT_ID SNPEFF_AMINO_ACID_CHANGE SNPEFF_CODON_CHANGE
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

SNP_FILTERS = --filterName 'QD' --filterExpression 'QD < 2.0' --filterName 'MQ' --filterExpression 'MQ < 40.0' --filterName 'FS' --filterExpression 'FS > 60.0' --filterName 'HapScore' --filterExpression 'HaplotypeScore > 13.0' --filterName 'MQRankSum' --filterExpression 'MQRankSum < -12.5' --filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < -8.0' --filterName 'Depth' --filterExpression 'DP < 5'
SNP_FILTER_INTRON = false
ifeq ($(SNP_FILTER_INTRON),true)
SNP_FILTERS += --filterName 'SnpEff' --filterExpression \"SNPEFF_EFFECT == 'UPSTREAM' || SNPEFF_EFFECT == 'DOWNSTREAM' || SNPEFF_EFFECT == 'INTERGENIC' || SNPEFF_EFFECT == 'INTRAGENIC' || SNPEFF_EFFECT == 'INTRON'\"
endif

TABLE_FIELDS = CHROM POS ID REF ALT AD DP PL AC AF AN BaseQRankSum DB Dels FS HaplotypeScore MLEAC MLEAF MQ MQ0 MQRankSum QD SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT SNPEFF_GENE_BIOTYPE SNPEFF_GENE_NAME SNPEFF_TRANSCRIPT_ID

GT_TABLE_FIELDS = GT AD DP PL

REPORT_STRATIFICATION := Filter
