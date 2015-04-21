ifndef CONFIG_INC
REF ?= hg19
SEQ_PLATFORM ?= illumina
SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt # generate using scripts/prepareFastq.sh [remove underscores from sample names]

##Set environment variables
export TMPDIR := $(HOME)/share/tmp
PATH := $(HOME)/share/usr/bin:$(PATH)
PATH := $(PATH):$(HOME)/share/usr/ensembl-git-tools/bin
export PATH
#PYTHONPATH := $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7:/opt/common/python/python-2.7.5/lib:/opt/common/python/python-2.7.5/lib/python2.7:$(HOME)/share/usr/anaconda/lib/python2.7/site-packages:$(HOME)/share/usr/anaconda/lib/python2.7:$(HOME)/share/usr/anaconda/lib/python
#export PYTHONPATH
export R_LIBS := $(HOME)/share/usr/lib64/R/library:/home/limr/R/x86_64-unknown-linux-gnu-library/3.0

LD_LIBRARY_PATH := $(HOME)/share/usr/mysql/lib:$(LD_LIBRARY_PATH)
export LD_LIBRARY_PATH

PERL5LIB := /opt/common/Perl/lib/perl5/:$(HOME)/share/usr/perl5/lib/perl5:$(HOME)/share/usr/ActivePerl-5.18/lib
PERL5LIB := $(HOME)/share/usr/src/bioperl-live:$(PERL5LIB)
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-compara/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-variation/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-funcgen/modules
export PERL5LIB

##Interpreters
PYTHON = /opt/common/python/python-2.7.5/bin/python
ANACONDA_PYTHON = $(HOME)/share/usr/anaconda/bin/python
PERL = $(HOME)/share/usr/ActivePerl-5.18/bin/perl
#RSCRIPT = /usr/bin/Rscript
#R = /opt/common/R/R-3.0.0/bin/R
R = $(HOME)/share/usr/bin/R-3.1.2
#RSCRIPT = /opt/common/R/R-3.0.0/bin/Rscript
MY_RSCRIPT = $(HOME)/share/usr/bin/Rscript-3.1.2
RSCRIPT = $(MY_RSCRIPT)
#RSCRIPT = $(HOME)/usr/bin/Rscript

JARDIR := $(HOME)/share/usr/lib/java

### Applications
MUTECT_JAR := $(JARDIR)/muTect-1.1.4.jar
#BOWTIE = /opt/common/bowtie2/bowtie2-2.1.0/bowtie2
BOWTIE = $(PERL) $(HOME)/share/usr/bowtie2-2.2.4/bowtie2
#BWA = /opt/common/bwa/bwa-0.6.2-r126/bwa
BWA = $(HOME)/share/usr/bin/bwa-0.7.10
#PBWA = mpiexec $(HOME)/share/usr/bin/pBWA
#GSNAP = $(HOME)/share/usr/bin/gsnap
#GSNAP_REF = $(HOME)/share/references/gmapdb/
#GSNAP_MAPS = $(GSNAP_REF)/$(REF)/$(REF).maps
SAMTOOLS = $(HOME)/share/usr/bin/samtools-0.1.19
SAMTOOLS2 = $(HOME)/share/usr/bin/samtools-1.1
VCFUTILS = $(HOME)/share/usr/bin/vcfutils.pl
BCFTOOLS = /opt/bin/bcftools
BCFTOOLS2 = $(HOME)/share/usr/bin/bcftools-1.1
BEDTOOLS = /opt/common/bedtools/bedtools-2.17.0/bin/bedtools
BGZIP = $(HOME)/share/usr/bin/bgzip
IGVTOOLS = $(HOME)/share/usr/IGVTools/igvtools
VCFTOOLS = $(HOME)/share/usr/bin/vcftools-0.1.10
VCF_SORT = $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl
SNP_EFF = $(JAVA) -Xmx4G -jar $(JARDIR)/snpEff.jar
CUFFLINKS = $(HOME)/share/usr/bin/cufflinks-2.2.1
CUFFCMP = $(HOME)/share/usr/bin/cuffcompare-2.2.1
#SNVMIX = $(HOME)/share/usr/bin/SNVMix2
#BAM2FASTQ = $(HOME)/share/usr/bin/bam2fastq
#QUALIMAP = $(HOME)/share/usr/qualimap/qualimap
TOPHAT = $(HOME)/share/usr/bin/tophat2-2.0.13
#DEFUSE = /opt/common/defuse/defuse-0.6.1/scripts/defuse.pl
DEFUSE = $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl

#VARSCAN
ONCOFUSE_JAR = $(JARDIR)/oncofuse-1.0.8.jar
VARSCAN_JAR = $(JARDIR)/VarScan.v2.3.7.jar

# PICARD tools
PICARD_DIR = $(JARDIR)/picard-1.92

# scripts
SCRIPTS_DIR = scripts
MERGE = $(SCRIPTS_DIR)/merge.R
VCF_TO_TABLE = $(SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP = $(SCRIPTS_DIR)/posnGeneLookup.pl
RBIND = $(SCRIPTS_DIR)/rbind.R
NORMAL_FILTER = $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl
SOMATIC_FILTER_VCF = $(PERL) $(SCRIPTS_DIR)/somaticFilterVCF.pl

#JAVA_TMP := $(HOME)/tmp
JAVA_BIN = /usr/lib/jvm/jre-1.6.0/bin/java
JAVA7_BIN = /usr/bin/java
#JAVA = /usr/bin/java $(JAVA_ARGS)

GET_INSERT_SIZE = $(HOME)/share/usr/bin/getInsertSize.py
endif
CONFIG_INC = true
