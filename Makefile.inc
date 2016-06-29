# vim: set ft=make:
ifneq ("$(wildcard config.inc)", "")
include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
include project_config.inc
endif
include modules/config.inc

# default reference: 
REF ?= hg19
include modules/$(REF).inc

ifndef MAKEFILE_INC

SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt # generate using scripts/prepareFastq.sh [remove underscores from sample names]

PROJECT_NAME = $(shell pwd | sed 's:.*/projects/::; s:.*/data/::; s:.*kinglab/::; s:/:_:g')

ifneq ($(wildcard $(SAMPLE_FILE)),)
  SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_FILE))
endif

get_tumors = $(wordlist 1,$(shell expr $(words $1) - 1),$1)
get_normal = $(lastword $1)
get_space_sets = $(shell sed '/^\#/d; s/\s\+/ /g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')
get_underscore_sets = $(shell sed '/^\#/d; s/\s\+/_/g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')

ifneq ($(wildcard $(SAMPLE_SET_FILE)),)
  NUM_SETS := $(shell sed '/^\#/d' $(SAMPLE_SET_FILE) | wc -l | cut -f 1 -d' ')
  SETS_SEQ := $(shell seq 1 $(NUM_SETS))
  $(foreach i,$(SETS_SEQ), \
	  $(eval set.$i := $(call get_space_sets,$i)))
  # define set_lookup.$(sample)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set_lookup.$(sample) := $i)))
  # define SAMPLE_SETS to contain sample_sample_normal
  $(foreach i,$(SETS_SEQ), \
	  $(eval SAMPLE_SETS += $(call get_underscore_sets,$i)))
  # tumor.SET => tumors
  $(foreach i,$(SETS_SEQ), \
	$(eval tumor.$(call get_underscore_sets,$i) := $(call get_tumors,$(set.$i))))
  # normal.SET => normal
  $(foreach i,$(SETS_SEQ), \
	$(eval normal.$(call get_underscore_sets,$i) := $(call get_normal,$(set.$i))))
  NORMAL_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_normal,$(set.$i)))
  TUMOR_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_tumors,$(set.$i)))
  SAMPLE_PAIRS = $(foreach set,$(SAMPLE_SETS),$(foreach tumor,$(tumor.$(set)),$(tumor)_$(normal.$(set))))
  # define tumor.pair and normal.pair to retrieve tumor/normal from pairs
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor.$(tumor)_$(normal.$(set)) := $(tumor))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval normal.$(tumor)_$(normal.$(set)) := $(normal.$(set)))))
  # tumor_normal.$(tumor) to retrieve tumor_normal
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(tumor) := $(tumor)_$(normal.$(set)))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(normal.$(set)) := $(tumor)_$(normal.$(set)))))
  # set.$(sample) to retrieve sample sets (underscores)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set.$(sample) := $(word $i,$(SAMPLE_SETS)))))
  UNMATCHED_SAMPLES = $(shell sed '/^\#/d' $(SAMPLE_FILE) $(SAMPLE_SET_FILE) | tr ' ' '\n' | sort | uniq -u)

  SAMPLE_SET_PAIRS = $(shell echo "$(SAMPLE_PAIRS) $(SAMPLE_SETS)" | tr ' ' '\n' | sort | uniq)
  $(foreach set,$(SAMPLE_SET_PAIRS), \
	$(eval samples.$(set) := $(shell echo "$(set)" | sed 's/_/ /g')))
endif

ifneq ($(wildcard $(SAMPLE_SPLIT_FILE)),)
  A = $(shell cut -f1 $(SAMPLE_SPLIT_FILE))
  B = $(shell cut -f2 $(SAMPLE_SPLIT_FILE))
  $(foreach i,$(shell seq 1 $(words $(A))),$(eval split.$(word $i,$(A)) += $(word $i,$(B))))
  UNSPLIT_SAMPLES = $(B)
  SPLIT_SAMPLES = $(shell cut -f1 $(SAMPLE_SPLIT_FILE) | sort | uniq)
endif


ALL_SAMPLES = $(SAMPLE_PAIRS) $(SAMPLES)

SHELL = /bin/bash

export TMPDIR := $(HOME)/share/tmp

#VARSCAN
VARSCAN_MEM = $(JAVA) -Xmx$1 -jar $(VARSCAN_JAR)
VARSCAN = $(call VARSCAN_MEM,8G)

# PICARD tools
PICARD_MEM = 10G
PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
PICARD_DIR ?= $(HOME)/share/usr/lib/java

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

PICARD = $(JAVA8) -Xmx$2 -jar $(PICARD2_JAR) $1 $(PICARD_OPTS)

# SNP EFF
SNP_EFF = $(call SNP_EFF_MEM,$(DEFAULT_JAVA_MEM))
SNP_EFF_MEM = $(JAVA7) -Xmx$1 -jar $(SNP_EFF_JAR)
SNP_SIFT_MEM = $(JAVA7) -Xmx$1 -jar $(SNP_SIFT_JAR)
SNP_SIFT = $(call SNP_SIFT_MEM,$(DEFAULT_JAVA_MEM))
VCF_EFF_ONE_PER_LINE = $(HOME)/share/usr/snpEff-4.1/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = modules/scripts/joinEff.pl

COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)


JAVA_ARGS = -Djava.io.tmpdir=$(TMPDIR)
#JAVA_TMP := $(HOME)/tmp
JAVA := $(JAVA_BIN) $(JAVA_ARGS)
JAVA6 := $(JAVA6_BIN) $(JAVA_ARGS)
JAVA7 := $(JAVA7_BIN) $(JAVA_ARGS)
JAVA8 := $(JAVA8_BIN) $(JAVA_ARGS)
#JAVA = /usr/bin/java $(JAVA_ARGS)

NO_RM ?= false
ifeq ($(NO_RM),true)
    RM := touch
    RMR = touch
else
    RM := rm -f
    RMR := rm -r
endif

# define $(,) and $( ) for escaping commas and spaces
, := ,
space :=
space +=
$(space) := 
$(space) +=



# $(call strip-suffix,filename)
strip-suffix = $(firstword $(subst ., ,$(1)))


LINK = ln -svf $(notdir $1) $(2) && touch $1

NOW := $(shell date +"%F")


MKDIR = mkdir -p -m 775
MKDIRS = $(MKDIR) $(LOGDIR)/$(@D) $(@D)
LOG = $(LOGDIR)/$(@).log



UMASK = 002


ACTIVATE_JRFLAB_MODULES_ENV = unset PYTHONPATH; source $(JRFLAB_MODULES_ENV)/bin/activate $(JRFLAB_MODULES_ENV)
UMASK = 002

# SGE variables
QUEUE ?= jrf.q all.q

### PBS
# mem for PBS: multiply mem_per_cpu by ncpu for request
# $(call PBS_MEM,ncpu,mem_per_cpu)
PBS_MEM = -l mem=$2
PBS_PARALLEL_MEM = -l nodes=1:ppn=$1 -l mem=$(shell unit=`echo $2 | sed 's/[0-9]//'`; echo "$1 * $2" | sed 's/[A-Za-z]//g' | bc | sed 's/$$/'$$unit'/')
PBS_QSUB_ARGS = -e $(LOG) -j oe

#### SGE
# mem for SGE
# $(call SGE_MEM,ncpu,mem_per_cpu)
SGE_PARALLEL_ENV = smp
SGE_PARALLEL_MEM = -pe $(SGE_PARALLEL_ENV) -l virtual_free=$1,h_vmem=$2
SGE_MEM = -l virtual_free=$1,h_vmem=$2
QSUB_PRIORITY ?= 0
SGE_QSUB_ARGS = -V -wd $(PWD) -now n -notify -o $(LOG) -j y -S $(SHELL) -p $(QSUB_PRIORITY) -b n


## DRMAA qsub script
QSUB_SCRIPT = perl modules/scripts/qsub.pl
QSUB_SCRIPT_OPTS = -o $@
QS = $(QSUB_SCRIPT) $(QSUB_SCRIPT_OPTS) $1 -- $(foreach q,$(strip $2), -q $q)
QSUB = $(call QS,,$(QUEUE))
QSUB_QUEUE = $(call QS,,$1)
# check output file for non-zero size
QSUB_CHECK = $(call QS,-c)
QSUB_CHECK_QUEUE = $(call QS,-c,-q $1)

#### no cluster
INIT = $(MKDIRS); umask $(UMASK); set -o pipefail; $(ACTIVATE_JRFLAB_MODULES_ENV);
INIT_NOCLUSTER = echo "umask $(UMASK); set -o pipefail; $(ACTIVATE_JRFLAB_MODULES_ENV);";

CLUSTER_ENGINE ?= SGE
USE_CLUSTER ?= true
ifeq ($(USE_CLUSTER),false)
LSCRIPT = $(MKDIRS); ($(INIT_NOCLUSTER) echo $1) | $(SHELL) >& $(LOG)
LSCRIPT_PARALLEL_MEM = $(MKDIRS); ($(INIT_NOCLUSTER) echo $4) | $(SHELL) >& $(LOG)
LSCRIPT_MEM = $(MKDIRS); ($(INIT_NOCLUSTER) echo $3) | $(SHELL) >& $(LOG)
LSCRIPT_NAMED = $(call LSCRIPT,$2)
LSCRIPT_NAMED_PARALLEL_MEM = $(call LSCRIPT_PARALLEL_MEM,$2,$3,$4,$5)
LSCRIPT_NAMED_MEM = $(call LSCRIPT_MEM,$2,$3,$4)
else
LSCRIPT = $(call LSCRIPT_NAMED_MEM,$(@F),1G,2G,$1)
LSCRIPT_PARALLEL_MEM = $(call LSCRIPT_NAMED_PARALLEL_MEM,$(@F),$1,$2,$3,$4)
LSCRIPT_MEM = $(call LSCRIPT_NAMED_MEM,$(@F),$1,$2,$3)

#$(call LS,qsub,name,qsub_args,script)
LS = $(MKDIRS); (echo "umask $(UMASK); set -o pipefail; $(ACTIVATE_JRFLAB_MODULES_ENV);" $4) | $1 -N $(PROJECT_NAME)_$2 $3
LSCRIPT_NAMED = $(call LS,$(QSUB),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,1G,2G),$2)
LSCRIPT_NAMED_PARALLEL_MEM = $(call LS,$(QSUB),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_PARALLEL_MEM,$2,$3,$4),$5)
LSCRIPT_NAMED_MEM = $(call LS,$(QSUB),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,$2,$3),$4)

LSCRIPT_CHECK = $(call LSCRIPT_CHECK_NAMED,$(@F),$1)
LSCRIPT_CHECK_PARALLEL_MEM = $(call LSCRIPT_CHECK_NAMED_PARALLEL_MEM,$(@F),$1,$2,$3,$4)
LSCRIPT_CHECK_MEM = $(call LSCRIPT_CHECK_NAMED_MEM,$(@F),$1,$2,$3)

LSCRIPT_CHECK_NAMED = $(call LS,$(QSUB_CHECK),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,1G,2G),$2)
LSCRIPT_CHECK_NAMED_PARALLEL_MEM = $(call LS,$(QSUB_CHECK),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_PARALLEL_MEM,$2,$3,$4),$5)
LSCRIPT_CHECK_NAMED_MEM = $(call LS,$(QSUB_CHECK),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,$2,$3),$4)
endif

#MD5 = md5sum $(1:.md5=) > $1
MD5 = md5sum $(@:.md5=) > $@
#CHECK_MD5 = md5sum -c $(filter %.md5,$^); i=0; until md5sum --status -c $(filter %.md5,$^); do if (( $$i % 5  == 0 )); then date; md5sum -c $(filter %.md5,$^); fi; if (( $$i > 20 )); then echo giving up; exit 1; fi; sleep 5; let i=i+1 ; done;
CHECK_MD5 = for i in {0..20}; do if md5sum -c $(filter %.md5,$^); then break; fi; sleep 5; done;
<M = $(<:.md5=)
^M = $(^:.md5=)
@M = $(@:.md5=)
<<M = $(word 2,$(^M))
<<<M = $(word 3,$(^M))

# $(call CHECK_VCF,vcf.file,target,command)
# check for valid header, ie non-empty vcf file
CHECK_VCF = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi
CHECK_MAF = if [ `grep -v '^\#' $1 | wc -l` -eq 1 ] && [ `grep '^Hugo_Symbol' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi

<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
4< = $(word 4,$^)
<<<<< = $(word 5,$^)
5< = $(word 5,$^)

KNIT = $(RSCRIPT) modules/scripts/knit.R


define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef

##Set environment variables
export TMPDIR := $(HOME)/share/tmp
PATH := $(HOME)/share/usr/bin:$(PATH)
PATH := $(PATH):$(HOME)/share/usr/ensembl-git-tools/bin
export PATH
export R_LIBS := $(HOME)/share/usr/lib64/R/library:/home/limr/R/x86_64-unknown-linux-gnu-library/3.0

LD_LIBRARY_PATH := $(HOME)/share/usr/mysql/lib:$(LD_LIBRARY_PATH)
export LD_LIBRARY_PATH

PERL5LIB := $(HOME)/share/usr/perl5/lib/perl5:$(HOME)/share/usr/ActivePerl-5.18/lib
PERL5LIB := $(HOME)/share/usr/src/bioperl-live:$(PERL5LIB)
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-compara/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-variation/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-funcgen/modules
export PERL5LIB

endif
MAKEFILE_INC = true
