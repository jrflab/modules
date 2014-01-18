# vim: set ft=make:

# default reference: 
REF ?= hg19
include $(HOME)/share/modules/$(REF).inc

#SAMPLE_PAIR_FILE ?= sample_pairs.txt
SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt # generate using scripts/prepareFastq.sh [remove underscores from sample names]

#ifneq ($(wildcard $(SAMPLE_PAIR_FILE)),)
#TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
#NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
#$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
#$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))
#SAMPLE_PAIRS ?= $(foreach tumor,$(TUMOR_SAMPLES),$(tumor)_$(normal_lookup.$(tumor)) )
#endif

ifneq ($(wildcard $(SAMPLE_FILE)),)
  SAMPLES ?= $(shell cat $(SAMPLE_FILE))
endif

get_tumors = $(wordlist 1,$(shell expr $(words $1) - 1),$1)
get_normal = $(lastword $1)

ifneq ($(wildcard $(SAMPLE_SET_FILE)),)
  NUM_SETS := $(shell wc -l $(SAMPLE_SET_FILE) | cut -f 1 -d' ')
  SETS_SEQ := $(shell seq 1 $(NUM_SETS))
  $(foreach i,$(SETS_SEQ),$(eval set.$i := $(shell sed -n $(i)p $(SAMPLE_SET_FILE))))
  $(foreach i,$(SETS_SEQ),$(foreach sample,$(shell sed -n $(i)p $(SAMPLE_SET_FILE)),$(eval set_lookup.$(sample) := $i)))
  ifndef SAMPLE_SETS
  $(foreach i,$(SETS_SEQ),$(eval SAMPLE_SETS += $(shell sed -n 's/ /_/g;$(i)p' $(SAMPLE_SET_FILE))))
  endif
  SAMPLE_PAIRS ?= $(foreach i,$(SETS_SEQ), \
				  $(foreach tumor,$(call get_tumors,$(set.$i)), \
				  $(tumor)_$(call get_normal,$(set.$i))))
  TUMOR_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_tumors,$(set.$i)))
  NORMAL_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_normal,$(set.$i)))

#ifeq ($(NORMAL_MATCHED),true)
#ifndef SAMPLE_PAIRS
#TUMOR_SAMPLES ?= $(shell sed 's/\w*$$//; s/ /\n/' $(SAMPLE_SET_FILE) | sort | uniq)
#NORMAL_SAMPLES ?= $(shell sed 's/.* //' $(SAMPLE_SET_FILE) | sort | uniq)
#$(foreach i,$(shell seq 1 $(NUM_SETS)),$(foreach tumor,$(wordlist 1,$(shell expr $(words $(set.$i)) - 1),$(set.$i)),$(eval normal_lookup.$(tumor) = $(lastword $(set.$i)))))
#$(foreach i,$(shell seq 1 $(NUM_SETS)),$(eval tumor_lookup.$(lastword $(set.$i)) += $(wordlist 1,$(shell expr $(words $(set.$i)) - 1),$(set.$i))))
#$(foreach tumor,$(TUMOR_SAMPLES),$(foreach normal,$(normal_lookup.$(tumor)),$(eval SAMPLE_PAIRS += $(tumor)_$(normal))))
#endif
endif




#$(foreach i,$(shell seq 1 $(NUM_SETS)),$(info $(wordlist 1,$(shell expr $(words $(set.$i)) - 1),$(set.$i))))


ifneq ($(wildcard $(SAMPLE_SPLIT_FILE)),)
A = $(shell cut -f1 $(SAMPLE_SPLIT_FILE))
B = $(shell cut -f2 $(SAMPLE_SPLIT_FILE))
$(foreach i,$(shell seq 1 $(words $(A))),$(eval split_lookup.$(word $i,$(A)) += $(word $i,$(B))))
SPLIT_SAMPLES = $(shell cut -f1 $(SAMPLE_SPLIT_FILE) | sort | uniq)
endif

#$(foreach i,$(shell seq 1 $(words $(LANES))),$(eval lane_lookup.$(word $i,$(LANES)) := $(word $i,$(LANE_SAMPLES))))
#$(foreach i,$(shell seq 1 $(words $(LANE_SAMPLES))),$(eval sample_lookup.$(word $i,$(LANE_SAMPLES)) += $(word $i,$(LANES))))


ALL_SAMPLES = $(SAMPLE_PAIRS) $(SAMPLES)

#SHELL = /bin/bash --noprofile -c 'umask 002; eval "$$1"'
SHELL = /bin/bash

#SHELL = umask $(UMASK); exec /bin/bash
#SHELL = export SGE_RREQ="-q all.q -now no" && umask $(UMASK) && exec /bin/bash

PATH := $(HOME)/share/usr/bin:$(PATH)
export PATH
PYTHONPATH := $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7:/opt/common/python/python-2.7.3/lib:/opt/common/python/python-2.7.3/lib/python2.7
export PYTHONPATH
#export PYTHONHOME := /opt/common/python/python-2.7.3/
#Set environment variables
#export R_LIBS_USER := $(HOME)/share/usr/lib64/R/library:${R_LIBS_USER}
export R_LIBS := $(HOME)/share/usr/lib64/R/library:/home/limr/R/x86_64-unknown-linux-gnu-library/3.0

LD_LIBRARY_PATH := $(HOME)/share/usr/mysql/lib:$(LD_LIBRARY_PATH)
export LD_LIBRARY_PATH

#PERL5LIB := $(HOME)/share/usr/lib/perl5/:${PERL5LIB}
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl/modules/
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl-variation/modules/
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl-functgenomics/modules/
#PERL5LIB += :$(HOME)/share/usr/lib/ensembl-compara/modules/

PERL5LIB := /opt/common/Perl/lib/perl5/

export PERL5LIB

#Interpreters
PYTHON = /opt/common/python/python-2.7.3/bin/python
PERL = /usr/bin/perl
#RSCRIPT = /usr/bin/Rscript
R = /opt/common/R/R-3.0.0/bin/R
RSCRIPT = /opt/common/R/R-3.0.0/bin/Rscript
#RSCRIPT = $(HOME)/usr/bin/Rscript

JARDIR := $(HOME)/share/usr/lib/java

### Applications
BOWTIE = /opt/common/bowtie2/bowtie2-2.1.0/bowtie2
#BWA = /opt/common/bwa/bwa-0.6.2-r126/bwa
BWA = $(HOME)/share/usr/bin/bwa
#PBWA = mpiexec $(HOME)/share/usr/bin/pBWA
#GSNAP = $(HOME)/share/usr/bin/gsnap
#GSNAP_REF = $(HOME)/share/references/gmapdb/
#GSNAP_MAPS = $(GSNAP_REF)/$(REF)/$(REF).maps
SAMTOOLS = $(HOME)/share/usr/bin/samtools
BCFTOOLS = /opt/bin/bcftools
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
DEFUSE = $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl
#DEFUSE = /opt/common/defuse/defuse-0.6.1/scripts/defuse.pl
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
NORMAL_FILTER = $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl

JAVA_ARGS = -Djava.io.tmpdir=$(TMPDIR)
#JAVA_TMP := $(HOME)/tmp
JAVA := /usr/lib/jvm/jre-1.6.0/bin/java $(JAVA_ARGS)
JAVA7 := /usr/bin/java $(JAVA_ARGS)
#JAVA = /usr/bin/java $(JAVA_ARGS)

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
CHR_CHROMOSOMES ?= false
ifdef INTERVALS_FILE
   CHROMOSOMES := $(shell grep -v '@' $(INTERVALS_FILE) | cut -f1 | sort | uniq)
else
   ifeq ($(CHR_CHROMOSOMES),true)
      CHROMOSOMES := chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrM
      ifeq ($(INCLUDE_CHR_Y),true)
	    CHROMOSOMES += chrY
      endif
   else
      CHROMOSOMES := 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X MT
      ifeq ($(INCLUDE_CHR_Y),true)
	    CHROMOSOMES += Y
      endif
   endif
endif
INCLUDE_CHR_Y ?= false




# $(call strip-suffix,filename)
strip-suffix = $(firstword $(subst ., ,$(1)))


LINK = ln -svf $(notdir $1) $(2) && touch $1

NOW := $(shell date +"%F")


MKDIR = mkdir -p -m 770
MKDIRS = $(MKDIR) $(LOGDIR)/$(@D) $(@D)
LOG = $(LOGDIR)/$(@:.md5=).log


#LOG_MEM = $(HOME)/share/scripts/monitorMem.sh $(LOG).top &
LOG_MEM = 

# SGE variables
QUEUE = jrf.q
SGE_RREQ = -q $(QUEUE) -N X$(@F) -V -cwd 
PARALLEL_ENV = smp
MEM_FREE = -l virtual_free=$1,h_vmem=$2
REQ = SGE_RREQ="$(SGE_RREQ)"

REQ_PARALLEL_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,$2,$3) -pe $(PARALLEL_ENV) $(1)"
REQ_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,$1,$2)"
UMASK = 002
INIT = $(MKDIRS); umask $(UMASK); $(LOG_MEM) 
# $(call INIT_MEM,mem_free,h_vmem)
INIT_MEM = $(MKDIRS); umask $(UMASK); $(call REQ_MEM,$1,$2); $(LOG_MEM) 
# $(call INIT_PARALLEL_MEM,ncpus,mem_free,h_vmem)
INIT_PARALLEL_MEM = $(MKDIRS); umask $(UMASK); $(call REQ_PARALLEL_MEM,$1,$2,$3); $(LOG_MEM) 

# openMPI
MPIRUN = /opt/bin/mpirun
MPI_ENV = openmpi
REQ_MPI_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,$2,$3) -pe $(MPI_ENV) $(1)"
INIT_MPI_MEM = $(MKDIRS); umask $(UMASK); $(call REQ_MPI_MEM,$1,$2,$3); $(MPIRUN) -np $1

QRSH = /common/sge/bin/lx-amd64/qrsh
LAUNCH = $(INIT) $(QRSH) -N X$(@F) -V -cwd -q $(QUEUE) -now n -notify
LAUNCH_PARALLEL_MEM = $(INIT) $(LAUNCH) -pe $(PARALLEL_ENV) $1 $(call MEM_FREE,$2,$3)
LAUNCH_MEM = $(INIT) $(LAUNCH) $(call MEM_FREE,$1,$2)

# launch read script from stdin
QSUB = /common/sge/bin/lx24-amd64/qsub
override QSUB_ARGS = -V -cwd -q $(QUEUE) -now n -sync y -notify -r y -e $(LOG).e -o $(LOG).o -S $(SHELL)
USE_CLUSTER ?= true
ifeq ($(USE_CLUSTER),false)
LSCRIPT = $(INIT) echo "umask $(UMASK);" $1 | $(SHELL) > $(LOG).o 2> $(LOG).e
LSCRIPT_PARALLEL_MEM = $(INIT) echo "umask $(UMASK);" $4 | $(SHELL) > $(LOG).o 2> $(LOG).e
LSCRIPT_MEM = $(INIT) echo "umask $(UMASK);" $3 | $(SHELL) > $(LOG).o 2> $(LOG).e
else
LSCRIPT = $(call LSCRIPT_NAMED,$(@F),$1)
LSCRIPT_PARALLEL_MEM = $(call LSCRIPT_NAMED_PARALLEL_MEM,$(@F),$1,$2,$3,$4)
LSCRIPT_MEM = $(call LSCRIPT_NAMED_MEM,$(@F),$1,$2,$3)

LSCRIPT_NAMED = $(INIT) echo "umask $(UMASK);" $2 | $(QSUB) -N X$1 $(QSUB_ARGS) -- -
LSCRIPT_NAMED_PARALLEL_MEM = $(INIT) echo "umask $(UMASK);" $5 | $(QSUB) -N X$1 $(QSUB_ARGS) -pe $(PARALLEL_ENV) $2 $(call MEM_FREE,$3,$4) -- -
LSCRIPT_NAMED_MEM = $(INIT) echo "umask $(UMASK);" $4 | $(QSUB) -N X$1 $(QSUB_ARGS) $(call MEM_FREE,$2,$3) -- -
endif

#MD5 = md5sum $(1:.md5=) > $1
MD5 = md5sum $(@:.md5=) > $@
#CHECK_MD5 = md5sum -c $(filter %.md5,$^); i=0; until md5sum --status -c $(filter %.md5,$^); do if (( $$i % 5  == 0 )); then date; md5sum -c $(filter %.md5,$^); fi; if (( $$i > 20 )); then echo giving up; exit 1; fi; sleep 5; let i=i+1 ; done;
CHECK_MD5 = for i in {0..20}; do if md5sum -c $(filter %.md5,$^); then break; fi; sleep 5; done;
<M = $(<:.md5=)
^M = $(^:.md5=)
@M = $(@:.md5=)

# $(call CHECK_VCF,vcf.file,target,command)
# check for valid header, ie non-empty vcf file
CHECK_VCF = if [ `grep -v '^\#' $1 | wc -l` -eq 0 && `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi

SOMATIC_FILTER_VCF = $(PERL) $(HOME)/share/scripts/somaticFilterVCF.pl

GET_INSERT_SIZE = $(HOME)/share/usr/bin/getInsertSize.py

#signature md5
#makepp_no_builtin = 1
#makepp_percent_subdirs := 1
