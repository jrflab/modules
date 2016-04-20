BAM_FIX_RG = false
BAM_NO_REALN = false
VCF_POST_ANN_FILTER_EXPRESSION = ExAC_AF > 0.05
ALIGNER = bwamem
BAM_NO_FILTER = false
BAM_DUP_TYPE = markdup
EXOME = true
GATK_HARD_FILTER_SNPS = true
BAM_REPROCESS = false
GATK_POOL_SNP_RECAL = false
BAM_NO_SORT = false
BAM_CHR1_BASE_RECAL = true
BAM_NO_RECAL = false
QSUB_PRIORITY = 0
BAM_PHRED64 = false
REF = b37

SAMPLES = GCT10T GCT11T GCT12T 
TUMOR_SAMPLES = GCT10T GCT11T GCT12T
NORMAL_SAMPLES = 
SAMPLE_PAIRS = 
SAMPLE_SETS = 

# tumor.normal = tumors

# normal.tumor = normal

# tumor.pair = tumor

# normal.pair = normal

# set.normal = set

# set.tumor = set

# set.pair = set

# sample_fastq_file
fq.GCT10T.0 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S7_L002_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S7_L002_R2_001.fastq.gz
fq.GCT10T.1 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S55_L005_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S55_L005_R2_001.fastq.gz
fq.GCT11T.0 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S54_L005_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S54_L005_R2_001.fastq.gz
fq.GCT11T.1 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S6_L002_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S6_L002_R2_001.fastq.gz
fq.GCT12T.0 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S56_L005_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S56_L005_R2_001.fastq.gz
fq.GCT12T.1 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S8_L002_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S8_L002_R2_001.fastq.gz
split.GCT10T.0 = GCT10T
split.GCT10T.1 = GCT10T
split.GCT10T = GCT10T.0 GCT10T.1 
split.GCT11T.0 = GCT11T
split.GCT11T.1 = GCT11T
split.GCT11T = GCT11T.0 GCT11T.1 
split.GCT12T.0 = GCT12T
split.GCT12T.1 = GCT12T
split.GCT12T = GCT12T.0 GCT12T.1 
SPLIT_SAMPLES = GCT11T.0 GCT11T.1 GCT12T.1 GCT12T.0 GCT10T.1 GCT10T.0

# defaults
GATK_JAR2 ?= $(JARDIR)/GenomeAnalysisTK-3.3-0.jar
JARDIR ?= $(HOME)/share/usr/lib/java
EMBL_MYSQLDB_PORT ?= 38493
FASTQC ?= $(PERL) $(HOME)/share/usr/FastQC/fastqc
MY_RSCRIPT ?= $(HOME)/share/usr/bin/Rscript-3.1.2
ANACONDA_PYTHON ?= $(HOME)/share/usr/anaconda/bin/python
VCFUTILS ?= $(HOME)/share/usr/bin/vcfutils.pl
UNZIP ?= /usr/bin/unzip
VT ?= $(HOME)/share/usr/bin/vt
BEDTOOLS ?= $(HOME)/share/usr/bin/bedtools
PERL ?= $(HOME)/share/usr/ActivePerl-5.18/bin/perl
VCF_SORT ?= $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl
TMPDIR ?= $(HOME)/share/tmp
SAMTOOLS2 ?= $(HOME)/share/usr/bin/samtools-1.2
RSCRIPT ?= $(HOME)/share/usr/bin/Rscript-3.1.2
GATK_JAR ?= $(JARDIR)/GenomeAnalysisTK.jar
SAMTOOLS ?= $(HOME)/share/usr/bin/samtools-0.1.19
SNP_EFF_JAR ?= $(JARDIR)/snpEff-4.1.jar
INTRON_POSN_LOOKUP ?= modules/scripts/posnGeneLookup.pl
JAVA8_BIN ?= $(HOME)/usr/jdk1.8.0_74/bin/java
EMBL_MYSQLDB_USER ?= embl
NSFP_FIELDS ?= Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_score Polyphen2_HVAR_pred 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF MutationAssessor_pred MutationAssessor_score MutationTaster_pred MutationTaster_score PROVEAN_pred ExAC_Adj_AF clinvar_rs clinvar_clnsig Interpro_domain
SOMATIC_FILTER_VCF ?= $(PERL) modules/scripts/somaticFilterVCF.pl
VARSCAN_JAR ?= $(JARDIR)/VarScan.v2.3.9.jar
BGZIP ?= $(HOME)/share/usr/bin/bgzip
MUTECT_JAR ?= $(JARDIR)/muTect-1.1.4.jar
INTEGRATE ?= $(HOME)/share/usr/Integrate_0.2.0/bin/Integrate
JAVA_BIN ?= /usr/lib/jvm/jre-1.6.0/bin/java
PICARD2_JAR ?= $(JARDIR)/picard-tools-2.1.1/picard.jar
SNP_EFF_CONFIG ?= $(HOME)/share/usr/snpEff-4.1/snpEff.config
EMBL_MYSQLDB_DB ?= homo_sapiens_core_75_37
TABIX ?= $(HOME)/share/usr/bin/tabix
ANACONDA_27_ENV ?= /home/debruiji/share/usr/anaconda-envs/anaconda-2.7
VCFTOOLS ?= $(HOME)/share/usr/bin/vcftools-0.1.10
EMBL_MYSQLDB_HOST ?= 10.0.200.71
JAVA6_BIN ?= /usr/lib/jvm/jre-1.6.0/bin/java
BOWTIE ?= $(PERL) $(HOME)/share/usr/bowtie2-2.2.4/bowtie2
BWA ?= $(HOME)/share/usr/bin/bwa-0.7.10
HISAT ?= $(HOME)/share/usr/hisat2-2.0.2-beta/hisat2
EMBL_MYSQLDB_PW ?= embl
BCFTOOLS2 ?= $(HOME)/share/usr/bin/bcftools-1.2
GET_INSERT_SIZE ?= $(HOME)/share/usr/bin/getInsertSize.py
SNP_SIFT_JAR ?= $(JARDIR)/snpSift-4.1.jar
CUFFLINKS ?= $(HOME)/share/usr/bin/cufflinks-2.2.1
DB_NSFP ?= $(HOME)/share/reference/snpEff-4.1/dbNSFP3.0b1a.hg19.txt.gz
JAVA7_BIN ?= /usr/bin/java
PICARD_DIR ?= $(JARDIR)/picard-1.92
PICARD_JAR ?= $(JARDIR)/picard-tools-1.141/picard.jar
DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl
CUFFCMP ?= $(HOME)/share/usr/bin/cuffcompare-2.2.1
VCF_TO_TABLE ?= modules/scripts/vcfToTable.R
NORMAL_FILTER ?= $(PERL) modules/scripts/normalFilterVCF.pl
MERGE ?= modules/scripts/merge.R
RBIND ?= modules/scripts/rbind.R
R ?= $(HOME)/share/usr/bin/R-3.1.2
SUFAM_ENV ?= /home/debruiji/share/usr/anaconda-envs/sufam-0.4.1
ONCOFUSE_JAR ?= $(HOME)/share/usr/oncofuse-1.0.9b2/Oncofuse.jar
PYTHON ?= $(HOME)/share/usr/anaconda/bin/python
TOPHAT ?= $(HOME)/share/usr/bin/tophat2-2.0.13
IGVTOOLS ?= $(HOME)/share/usr/IGVTools/igvtools
BCFTOOLS ?= $(HOME)/share/usr/bin/bcftools-0.1.19
