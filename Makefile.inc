# vim: set ft=make:

# default reference: 
REF ?= hg19

ifeq ($(REF),hg19)
	include ~/share/modules/hg19.inc
else ifeq ($(REF),mm10)
	include ~/share/modules/mm10.inc
endif

UMASK = 002
SHELL = /bin/bash --noprofile -c 'umask 002; eval "$$1"'

#SHELL = umask $(UMASK); exec /bin/bash
#SHELL = export SGE_RREQ="-q all.q -now no" && umask $(UMASK) && exec /bin/bash

export PATH := $(HOME)/share/usr/bin:$(PATH)
export PYTHONPATH := $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7:/opt/common/python/python-2.7.3/lib:/opt/common/python/python-2.7.3/lib/python2.7:$(PYTHONPATH)
export PYTHONHOME := /opt/common/python/python-2.7.3/
#Set environment variables
#export R_LIBS_USER := $(HOME)/share/usr/lib64/R/library:${R_LIBS_USER}
export R_LIBS := $(HOME)/share/usr/lib64/R/library

#PERL5LIB := $(HOME)/share/usr/lib/perl5/:${PERL5LIB}
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl/modules/
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl-variation/modules/
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl-functgenomics/modules/
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl-compara/modules/

export PERL5LIB

#Interpreters
PYTHON = /opt/common/python/python-2.7.3/bin/python
PERL = /usr/bin/perl
RSCRIPT = /usr/bin/Rscript
#RSCRIPT = $(HOME)/usr/bin/Rscript

JARDIR := $(HOME)/share/usr/lib/java

### Applications
BOWTIE = /opt/common/bowtie2/bowtie2-2.1.0/bowtie2
BWA = /opt/common/bwa/bwa-0.6.2-r126/bwa
#PBWA = mpiexec $(HOME)/share/usr/bin/pBWA
#GSNAP = $(HOME)/share/usr/bin/gsnap
#GSNAP_REF = $(HOME)/share/references/gmapdb/
#GSNAP_MAPS = $(GSNAP_REF)/$(REF)/$(REF).maps
SAMTOOLS = $(HOME)/share/usr/bin/samtools
IGVTOOLS = $(HOME)/share/usr/bin/IGVTools/igvtools
VCFTOOLS = $(HOME)/share/usr/bin/vcftools
VCF_SORT = $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl
SNP_EFF = $(JAVA) -Xmx4G -jar $(JARDIR)/snpEff.jar
SNP_EFF_CONFIG = $(HOME)/share/modules/snpEff.config 
CUFFLINKS = $(HOME)/share/usr/bin/cufflinks
CUFFCMP = $(HOME)/share/usr/bin/cuffcompare
#SNVMIX = $(HOME)/share/usr/bin/SNVMix2
#BAM2FASTQ = $(HOME)/share/usr/bin/bam2fastq
#QUALIMAP = $(HOME)/share/usr/qualimap/qualimap
#TOPHAT = $(HOME)/share/usr/tophat-2.0.8.Linux_x86_64/tophat2
#DEFUSE = $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl
TMPDIR = $(HOME)/share/tmp

# PICARD tools
PICARD_MEM = 10G
PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
PICARD_DIR = $(HOME)/share/usr/lib/java

ANALYZE_COV = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/AnalyzeCovariates.jar $(PICARD_OPTS)
SORT_SAM = $(call SORT_SAM_MEM,$(PICARD_MEM))
SORT_SAM_MEM = $(JAVA) -Xmx$1 -jar $(PICARD_DIR)/SortSam.jar $(PICARD_OPTS) TMP_DIR=$(TMPDIR)
REORDER_SAM = $(call REORDER_SAM_MEM,$(PICARD_MEM))
REORDER_SAM_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/ReorderSam.jar $(PICARD_OPTS)
MARK_DUP = $(call $(MARK_DUP_MEM,$(PICARD_MEM)))
MARK_DUP_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/MarkDuplicates.jar $(PICARD_OPTS)
MERGE_SAMS = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/MergeSamFiles.jar $(PICARD_OPTS)
INTERVAL_LIST_TOOL = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/IntervalListTools.jar $(PICARD_OPTS)
CREATE_SEQ_DICT = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CreateSequenceDictionary.jar $(PICARD_OPTS)
CALC_HS_METRICS = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CalculateHsMetrics.jar $(PICARD_OPTS)
COLLECT_MULT_METRICS = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CollectMultipleMetrics.jar $(PICARD_OPTS)
COLLECT_TARGETED_METRICS = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CollectTargetedPcrMetrics.jar $(PICARD_OPTS)

FIX_MATE = $(call FIX_MATE_MEM,$(PICARD_MEM))
FIX_MATE_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/FixMateInformation.jar $(PICARD_OPTS) TMP_DIR=$(TMPDIR)
SAM_TO_FASTQ = $(call SAM_TO_FASTQ_MEM,$(PICARD_MEM))
SAM_TO_FASTQ_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/SamToFastq.jar $(PICARD_OPTS)
CLEANBAM = $(call CLEANBAM_MEM,$(PICARD_MEM))
CLEANBAM_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/CleanSam.jar $(PICARD_OPTS)

ADD_RG = $(call ADD_RG_MEM,$(PICARD_MEM))
ADD_RG_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/AddOrReplaceReadGroups.jar $(PICARD_OPTS)
VALIDATE_SAM = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/ValidateSamFile.jar $(PICARD_OPTS)

# scripts
SCRIPTS_DIR = $(HOME)/share/scripts
MERGE = $(SCRIPTS_DIR)/merge.R
VCF_TO_TABLE = $(SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP = $(SCRIPTS_DIR)/posnGeneLookup.pl
RBIND = $(SCRIPTS_DIR)/rbind.R

JAVA_ARGS = -D.java.io.tmpdir=$(TMPDIR)
#JAVA_TMP := $(HOME)/tmp
#JAVA := /usr/lib/jvm/jre-1.6.0/bin/java $(JAVA_ARGS)
JAVA = /usr/bin/java $(JAVA_ARGS)

NO_RM ?= false
ifeq ($(NO_RM),true)
    RM := touch
else
    RM := rm -f
endif

# define $(,) and $( ) for escaping commas and spaces
, := ,
space :=
space +=
$(space) := 
$(space) +=

#CHROMOSOMES = chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrM chrX chrY
CHROMOSOMES = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT

# $(call strip-suffix,filename)
strip-suffix = $(firstword $(subst ., ,$(1)))


LINK = ln -svf $(notdir $1) $(2) && touch $1

NOW := $(shell date +"%F")


MKDIR = mkdir -p -m 770
MKDIRS = $(MKDIR) $(LOGDIR)/$(@D) $(@D)
LOG = $(LOGDIR)/$@.log

PARALLEL_ENV = smp
MEM_FREE = -l virtual_free=$1,h_vmem=$2

# SGE variables
SGE_RREQ = -q jrf.q -N X$(@F) -V -cwd
REQ = SGE_RREQ="$(SGE_RREQ)"
REQ_PARALLEL_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,$2,$3) -pe $(PARALLEL_ENV) $(1)"
REQ_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,$1,$2)"
INIT = $(REQ) $(MKDIRS);
# $(call INIT_MEM,mem_free,h_vmem)
INIT_MEM = $(MKDIRS); $(call REQ_MEM,$1,$2)
# $(call INIT_PARALLEL_MEM,ncpus,mem_free,h_vmem)
INIT_PARALLEL_MEM = $(MKDIRS); $(call REQ_PARALLEL_MEM,$1,$2,$3)

