ifndef CONFIG_INC

REF ?= b37
SEQ_PLATFORM ?= illumina
SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt # generate using scripts/prepareFastq.sh [remove underscores from sample names]

##Interpreters
ANACONDA_PYTHON ?= python
PYTHON ?= python
PERL ?= perl
R ?= R
MY_RSCRIPT ?= Rscript
RSCRIPT ?= Rscript

ANACONDA_27_ENV = $(HOME)/share/usr/anaconda-envs/anaconda-2.7
SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
MUTSIG_REPORT_ENV = $(HOME)/share/usr/anaconda-envs/mutsig-report-0.0.1
JRFLAB_MODULES_ENV = $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.4
ONCOTATOR_ENV = $(HOME)/share/usr/venv/oncotator-1.9.2.0
VEP_ENV = $(HOME)/share/usr/anaconda-envs/variant-effect-predictor-86
ASCAT_ENV = $(HOME)/share/usr/anaconda-envs/ascat
INNOVATION_ENV = $(HOME)/share/usr/env/innovation-lab-0.0.1
PIGZ_ENV ?= $(HOME)/share/usr/env/pigz-2.6
KALLISTO_ENV ?= $(HOME)/share/usr/env/kallisto-0.46.2
IMMUNE_ENV ?= $(HOME)/share/usr/env/r-immunedeconv-2.1.0
SUMREADS_ENV = $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.6
PYCLONE_ENV = $(HOME)/share/usr/env/pyclone-vi-0.1.2
PYCLONE_13_ENV = $(HOME)/share/usr/env/pyclone-0.13.1
GRIDSS_ENV = $(HOME)/share/usr/env/gridss-2.13.2
SVABA_ENV ?= $(HOME)/share/usr/env/svaba-1.1.0
SURVIVOR_ENV ?= $(HOME)/share/usr/env/survivor-1.0.7
ANNOTATE_SV_ENV ?= $(HOME)/share/usr/env/annot_sv-3.1.3
VIOLA_ENV = $(HOME)/share/usr/env/viola-sv-1.0.2
SIGNATURE_TOOLS_ENV = $(HOME)/share/usr/env/r-signature.tools.lib-2.2.0
CNVKIT_ENV ?= $(HOME)/share/usr/env/cnvkit-0.9.9
STARFISH_ENV ?= $(HOME)/share/usr/env/r-starfish-0.11
MEDICC_ENV = $(HOME)/share/usr/env/medicc2-0.8.1
VARIANT_ANNOTATION_ENV = $(HOME)/share/usr/env/r-variantannotation-1.44.0
FACETS_SUITE_ENV = $(HOME)/share/usr/env/r-facets-suite-2.0.8
CRAVAT_ENV = $(HOME)/share/usr/anaconda-envs/open-cravat
POLYSOLVER_ENV = $(HOME)/share/usr/anaconda-envs/hla-polysolver

JARDIR ?= $(HOME)/share/usr/lib/java

## Applications
UNZIP ?= /usr/bin/unzip
FASTQC ?= $(PERL) $(HOME)/share/usr/FastQC/fastqc
MUTECT_JAR ?= $(JARDIR)/muTect-1.1.7.jar
BOWTIE ?= $(PERL) $(HOME)/share/usr/bowtie2-2.2.4/bowtie2
BWA ?= bwa
SAMTOOLS ?= samtools
SAMTOOLS2 ?= samtools
VCFUTILS ?= $(HOME)/share/usr/bin/vcfutils.pl
BCFTOOLS2 ?= bcftools
BCFTOOLS ?= bcftools
PIGZ ?= pigz
BEDTOOLS ?= $(HOME)/share/usr/bin/bedtools
BGZIP ?= $(HOME)/share/usr/bin/bgzip
IGVTOOLS ?= $(HOME)/share/usr/IGVTools/igvtools
VCFTOOLS ?= $(HOME)/share/usr/bin/vcftools-0.1.10
VCF_SORT ?= $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl
SNP_EFF_JAR ?= $(JARDIR)/snpEff-4.3.jar
SNP_SIFT_JAR ?= $(JARDIR)/SnpSift-4.3.jar
SNP_EFF_CONFIG ?= modules/config/snpEff.conf
DB_NSFP ?= $(HOME)/share/reference/snpEff-4.1/dbNSFP3.0b1a.hg19.txt.gz
NSFP_FIELDS ?= Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_score Polyphen2_HVAR_pred 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF MutationAssessor_pred MutationAssessor_score MutationTaster_pred MutationTaster_score PROVEAN_pred ExAC_Adj_AF clinvar_rs clinvar_clnsig Interpro_domain
CUFFLINKS ?= cufflinks
CUFFCMP ?= cuffcompare
TOPHAT ?= tophat
DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl
ONCOFUSE_JAR ?= $(HOME)/share/usr/oncofuse-1.0.9b2/Oncofuse.jar
VARSCAN_JAR ?= $(JARDIR)/VarScan.v2.3.9.jar
MEDICC ?= medicc2

## PICARD tools
PICARD_DIR ?= $(JARDIR)/picard-1.92
PICARD_JAR ?= $(JARDIR)/picard-tools-1.141/picard.jar

## scripts
SCRIPTS_DIR ?= modules/scripts
MERGE ?= $(SCRIPTS_DIR)/merge.R
VCF_TO_TABLE ?= $(SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP ?= $(SCRIPTS_DIR)/posnGeneLookup.pl
RBIND ?= $(SCRIPTS_DIR)/rbind.R
NORMAL_FILTER ?= $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl
SOMATIC_FILTER_VCF ?= $(PERL) $(SCRIPTS_DIR)/somaticFilterVCF.pl
SUM_READS_RSCRIPT = $(RSCRIPT) $(SCRIPTS_DIR)/summarize_rnaseqreads.R
SUM_EXONS_RSCRIPT = $(RSCRIPT) $(SCRIPTS_DIR)/summarize_rnaseqreads_byexon.R
SUM_INTRONS_RSCRIPT = $(RSCRIPT) $(SCRIPTS_DIR)/summarize_rnaseqreads_byintron.R

JAVA_BIN ?=  $(JAVA8_BIN)
JAVA6_BIN ?= $(HOME)/share/usr/jdk1.6.0_45/bin/java
JAVA7_BIN ?= $(HOME)/share/usr/jdk1.7.0_45/bin/java
JAVA8_BIN ?= $(HOME)/share/usr/jdk1.8.0_121/bin/java

GET_INSERT_SIZE ?= $(HOME)/share/usr/bin/getInsertSize.py

## GATK
GATK_JAR ?= $(JARDIR)/GenomeAnalysisTK.jar
GATK_JAR2 ?= $(JARDIR)/GenomeAnalysisTK-3.7.jar

VT ?= $(HOME)/share/usr/bin/vt
TABIX ?= $(HOME)/share/usr/bin/tabix

CLUSTER_ENGINE ?= SGE
ifeq ($(CLUSTER_ENGINE),PBS)
MYSQLDB_HOST = cpu-6-2
EMBL_MYSQLDB_HOST = cpu-6-2
EMBL_MYSQLDB_PORT = 9992
EMBL_MYSQLDB_USER = embl
EMBL_MYSQLDB_PW = embl
EMBL_MYSQLDB_DB = homo_sapiens_core_85_37
endif
ifeq ($(CLUSTER_ENGINE),LSF)
MYSQLDB_HOST = swan
EMBL_MYSQLDB_HOST = swan
EMBL_MYSQLDB_PORT = 38493
EMBL_MYSQLDB_USER = embl
EMBL_MYSQLDB_PW = embl
EMBL_MYSQLDB_DB = homo_sapiens_core_75_37
endif
ifeq ($(CLUSTER_ENGINE),SGE)
MYSQLDB_HOST = ika
EMBL_MYSQLDB_HOST = ika
EMBL_MYSQLDB_PORT = 38493
EMBL_MYSQLDB_USER = embl
EMBL_MYSQLDB_PW = embl
EMBL_MYSQLDB_DB = homo_sapiens_core_75_37
endif
SNVBOX_CONF = $(PWD)/modules/db/snv_box-$(MYSQLDB_HOST).conf
FATHMM_CONF = $(PWD)/modules/db/fathmm_config-$(MYSQLDB_HOST).ini

INTEGRATE ?= $(HOME)/share/usr/Integrate_0.2.0/bin/Integrate

VCF2MAF ?= $(PERL) $(HOME)/share/usr/bin/vcf2maf.pl
MAF2MAF ?= $(PERL) $(HOME)/share/usr/bin/maf2maf.pl
MAF2VCF ?= $(PERL) $(HOME)/share/usr/bin/maf2vcf.pl

HISAT ?= $(HOME)/share/usr/hisat2-2.0.2-beta/hisat2

TVC ?= $(HOME)/share/usr/bin/tvc

ANNOVAR = $(PERL) $(HOME)/share/usr/annovar-2017-07-16/table_annovar.pl

VEP_PATH = $(VEP_ENV)/bin

SPLIT_BED = python modules/scripts/split_bed.py

INDEL_FILTER_VCF = python modules/vcf_tools/indel_filter_vcf.py
SNP_FILTER_VCF = python modules/vcf_tools/snp_filter_vcf.py

MERGE_VCF = python modules/vcf_tools/merge_vcf.py

# gbc command line
GBC ?= $(HOME)/share/usr/GetBaseCounts/GetBaseCounts

endif
CONFIG_INC = true
