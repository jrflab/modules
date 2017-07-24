# vim: set ft=make:
ifneq ("$(wildcard config.inc)", "")
include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
include project_config.inc
endif
include modules/config.inc
ifdef PROJECT_CONFIG
include $(PROJECT_CONFIG)
endif

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

VARSCAN = $(JAVA8) -Xmx10G -jar $(VARSCAN_JAR)

# SNP EFF
SNP_EFF = $(call SNP_EFF_MEM,$(DEFAULT_JAVA_MEM))
SNP_EFF_MEM = $(JAVA8) -Xmx$1 -jar $(SNP_EFF_JAR)
SNP_SIFT_MEM = $(JAVA8) -Xmx$1 -jar $(SNP_SIFT_JAR)
SNP_SIFT = $(call SNP_SIFT_MEM,$(DEFAULT_JAVA_MEM))
VCF_EFF_ONE_PER_LINE = $(HOME)/share/usr/snpEff-4.1/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = modules/scripts/joinEff.pl

COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)


JAVA_ARGS = -Djava.io.tmpdir=$(TMPDIR)
#JAVA_TMP := $(HOME)/tmp
JAVA := java $(JAVA_ARGS)
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
LOG = $(PWD)/$(LOGDIR)/$(@).log

UMASK = 002

export PATH := $(JRFLAB_MODULES_ENV)/bin:$(PATH)
ECHO_ACTIVATE_ENV = unset PYTHONPATH; OLDPATH=\$$PATH; tries=0; until [[ \$$tries -gt 10 ]] || source $1/bin/activate $1; do export PATH=\$$OLDPATH; tries=\`expr \$$tries + 1\`; jitter=\`expr \$$RANDOM % 10 + 1\`; sleep \`expr \$$jitter \\* \$$tries\`; done
ACTIVATE_ENV = unset PYTHONPATH; OLDPATH=$$PATH; tries=0; until [[ $$tries -gt 10 ]] || source $1/bin/activate $1; do export PATH=$$OLDPATH; tries=`expr $$tries + 1`; jitter=`expr $$RANDOM % 10 + 1`; sleep `expr $$jitter \* $$tries`; done

# SGE variables
#QUEUE ?= jrf.q all.q
QUEUE ?=

WALLTIMES = MEM SHORT LONG
### PBS
# mem for PBS: multiply mem_per_cpu by ncpu for request
# $(call PBS_MEM,ncpu,mem_per_cpu)
PBS_SHORT = -l mem=$2 -l walltime=1:00:00
PBS_MEM = -l mem=$2 -l walltime=24:00:00
PBS_LONG = -l mem=$2 -l walltime=128:00:00
PBS_PARALLEL_MEM_WALLTIME = -l nodes=1:ppn=$1 -l walltime=$3  -l mem=$(shell unit=`echo $2 | sed 's/[^A-Za-z]//g'`; echo "$1 * $2" | sed 's/[A-Za-z]//g; s/\..*//;' | bc | sed 's/$$/'$$unit'/')
PBS_PARALLEL_MEM = $(call PBS_PARALLEL_MEM_WALLTIME,$1,$3,12:00:00)
PBS_PARALLEL_SHORT = $(call PBS_PARALLEL_MEM_WALLTIME,$1,$3,1:00:00)
PBS_PARALLEL_LONG = $(call PBS_PARALLEL_MEM_WALLTIME,$1,$3,256:00:00)
PBS_QSUB_ARGS = -o $(LOG) -j oe -d $(PWD)
PBS_QSUB_SCRIPT = python modules/scripts/qsub_pbs.py
PBS_DOCKER_MEM = -l mem=$2 -l walltime=24:00:00 -l nodes=1:docker -l mem=$2
PBS_DOCKER_SHORT = -l mem=$2 -l walltime=1:00:00 -l nodes=1:docker -l mem=$2
PBS_DOCKER_LONG = -l mem=$2 -l walltime=128:00:00 -l nodes=1:docker -l mem=$2
PBS_DOCKER_PARALLEL_MEM_WALLTIME = -l nodes=1:docker:ppn=$1 -l walltime=$3 -l mem=$(shell unit=`echo $2 | sed 's/^[A-Za-z]//g'`; echo "$1 * $2" | sed 's/[A-Za-z]//g; s/\..*//;' | bc | sed 's/$$/'$$unit'/')
PBS_DOCKER_PARALLEL_MEM = $(call PBS_DOCKER_PARALLEL_MEM_WALLTIME,$1,$3,12:00:00)
PBS_DOCKER_PARALLEL_SHORT = $(call PBS_DOCKER_PARALLEL_MEM_WALLTIME,$1,$3,1:00:00)
PBS_DOCKER_PARALLEL_LONG = $(call PBS_DOCKER_PARALLEL_MEM_WALLTIME,$1,$3,128:00:00)
PBS_NODES = gpu-2-14 cpu-6-1 gpu-1-4


#### SGE
# mem for SGE
# $(call SGE_MEM,ncpu,mem_per_cpu)
SGE_PARALLEL_ENV = smp
SGE_PARALLEL_MEM = -pe $(SGE_PARALLEL_ENV) $1 -l virtual_free=$2,h_vmem=$3
SGE_PARALLEL_SHORT = $(SGE_PARALLEL_MEM)
SGE_PARALLEL_LONG = $(SGE_PARALLEL_MEM)
SGE_SHORT = $(SGE_MEM)
SGE_MEM = -l virtual_free=$1,h_vmem=$2
SGE_LONG = $(SGE_MEM)
QSUB_PRIORITY ?= 0
SGE_QSUB_ARGS = -V -wd $(PWD) -now n -notify -o $(LOG) -j y -S $(SHELL) -p $(QSUB_PRIORITY) -b n
SGE_NODES = e01 e02 e06
SGE_QSUB_SCRIPT = python modules/scripts/qsub.py



## DRMAA qsub script
QSUB_SCRIPT_OPTS = -o $@ -s $($(CLUSTER_ENGINE)_NODES)
QS = $($(CLUSTER_ENGINE)_QSUB_SCRIPT) $(QSUB_SCRIPT_OPTS) $1 --args $(foreach q,$(strip $2), -q $q)
QSUB = $(call QS,,$(QUEUE))
QSUB_QUEUE = $(call QS,,$1)
# check output file for non-zero size
QSUB_CHECK = $(call QS,-c)
QSUB_CHECK_QUEUE = $(call QS,-c,-q $1)


MAKE_PID := $(shell echo $$PPID)
JOB_FLAG := $(filter -j%, $(subst -j ,-j,$(shell ps T | grep "^\s*$(MAKE_PID).*$(MAKE)")))
JOBS := $(subst -j,,$(JOB_FLAG))

#### no cluster
RUN_FIRST = umask $(UMASK); $(if $(and $(findstring true,$(USE_CLUSTER)),$(JOBS)),sleep $$[ ( $$RANDOM % 20 ) + 1 ]s;) $(MKDIRS)
ECHO_FIRST = echo -n "umask $(UMASK); set -o pipefail; $(call ECHO_ACTIVATE_ENV,$(JRFLAB_MODULES_ENV)); ";
ECHO_FIRST_ENV = echo -n "umask $(UMASK); set -o pipefail; $(call ECHO_ACTIVATE_ENV,$1); ";

# run locally
INIT = $(RUN_FIRST); $(call ACTIVATE_ENV,$(JRFLAB_MODULES_ENV)); umask $(UMASK); set -o pipefail;
INIT_ENV = $(RUN_FIRST); $(call ACTIVATE_ENV,$1); umask $(UMASK); set -o pipefail;

CLUSTER_ENGINE ?= SGE
USE_CLUSTER ?= true
ifeq ($(USE_CLUSTER),false)
LSCRIPT = $(RUN_FIRST); ($(ECHO_FIRST) echo $1) | $(SHELL) >& $(LOG)
LSCRIPT_ENV = $(RUN_FIRST); ($(call ECHO_FIRST_ENV,$1)) echo $2) | $(SHELL) >& $(LOG)
LSCRIPT_PARALLEL_MEM = $(RUN_FIRST); ($(ECHO_FIRST) echo $4) | $(SHELL) >& $(LOG)
LSCRIPT_ENV_PARALLEL_MEM = $(RUN_FIRST); ($(call ECHO_FIRST_ENV,$1) echo $5) | $(SHELL) >& $(LOG)
LSCRIPT_PARALLEL_SHORT = $(LSCRIPT_PARALLEL_MEM)
LSCRIPT_ENV_PARALLEL_SHORT = $(LSCRIPT_ENV_PARALLEL_MEM)
LSCRIPT_PARALLEL_LONG = $(LSCRIPT_PARALLEL_MEM)
LSCRIPT_ENV_PARALLEL_LONG = $(LSCRIPT_ENV_PARALLEL_MEM)
LSCRIPT_MEM = $(RUN_FIRST); ($(ECHO_FIRST) echo $3) | $(SHELL) >& $(LOG)
LSCRIPT_ENV_MEM = $(RUN_FIRST); ($(call ECHO_FIRST_ENV,$1) echo $4) | $(SHELL) >& $(LOG)
LSCRIPT_SHORT = $(LSCRIPT_MEM)
LSCRIPT_ENV_SHORT = $(LSCRIPT_ENV_MEM)
LSCRIPT_LONG = $(LSCRIPT_MEM)
LSCRIPT_ENV_LONG = $(LSCRIPT_ENV_MEM)
LSCRIPT_CHECK = $(RUN_FIRST); ($(ECHO_FIRST) echo $1) | $(SHELL) >& $(LOG)
LSCRIPT_ENV_CHECK = $(RUN_FIRST); ($(call ECHO_FIRST_ENV,$1) echo $2) | $(SHELL) >& $(LOG)
LSCRIPT_CHECK_PARALLEL_MEM = $(RUN_FIRST) && ($(ECHO_FIRST) echo $4) | $(SHELL) >& $(LOG)
LSCRIPT_ENV_CHECK_PARALLEL_MEM = $(RUN_FIRST) && ($(call ECHO_FIRST_ENV,$1) echo $5) | $(SHELL) >& $(LOG)
LSCRIPT_CHECK_PARALLEL_SHORT = $(LSCRIPT_CHECK_PARALLEL_MEM)
LSCRIPT_ENV_CHECK_PARALLEL_SHORT = $(LSCRIPT_ENV_CHECK_PARALLEL_MEM)
LSCRIPT_CHECK_PARALLEL_LONG = $(LSCRIPT_CHECK_PARALLEL_MEM)
LSCRIPT_ENV_CHECK_PARALLEL_LONG = $(LSCRIPT_ENV_CHECK_PARALLEL_MEM)
LSCRIPT_CHECK_MEM = $(RUN_FIRST); ($(ECHO_FIRST) echo $3) | $(SHELL) >& $(LOG)
LSCRIPT_ENV_CHECK_MEM = $(RUN_FIRST); ($(call ECHO_FIRST_ENV,$1) echo $4) | $(SHELL) >& $(LOG)
LSCRIPT_CHECK_SHORT = $(LSCRIPT_CHECK_MEM)
LSCRIPT_ENV_CHECK_SHORT = $(LSCRIPT_ENV_CHECK_MEM)
LSCRIPT_CHECK_LONG = $(LSCRIPT_CHECK_MEM)
LSCRIPT_ENV_CHECK_LONG = $(LSCRIPT_ENV_CHECK_MEM)
LSCRIPT_NAMED = $(call LSCRIPT,$2)
LSCRIPT_ENV_NAMED = $(call LSCRIPT_ENV,$1,$3)
LSCRIPT_NAMED_PARALLEL_MEM = $(call LSCRIPT_PARALLEL_MEM,$2,$3,$4,$5)
LSCRIPT_ENV_NAMED_PARALLEL_MEM = $(call LSCRIPT_ENV_PARALLEL_MEM,$1,$3,$4,$5,$6)
LSCRIPT_NAMED_PARALLEL_SHORT = $(call LSCRIPT_PARALLEL_SHORT,$2,$3,$4,$5)
LSCRIPT_ENV_NAMED_PARALLEL_SHORT = $(call LSCRIPT_ENV_PARALLEL_MEM,$1,$3,$4,$5,$6)
LSCRIPT_NAMED_PARALLEL_LONG = $(call LSCRIPT_PARALLEL_LONG,$2,$3,$4,$5)
LSCRIPT_ENV_NAMED_PARALLEL_LONG = $(call LSCRIPT_ENV_PARALLEL_MEM,$1,$3,$4,$5,$6)
LSCRIPT_NAMED_MEM = $(call LSCRIPT_MEM,$2,$3,$4)
LSCRIPT_ENV_NAMED_MEM = $(call LSCRIPT_ENV_MEM,$1,$3,$4,$5)
LSCRIPT_NAMED_SHORT = $(LSCRIPT_NAMED_MEM)
LSCRIPT_ENV_NAMED_SHORT = $(LSCRIPT_ENV_NAMED_MEM)
LSCRIPT_NAMED_LONG = $(LSCRIPT_NAMED_MEM)
LSCRIPT_ENV_NAMED_LONG = $(LSCRIPT_ENV_NAMED_MEM)
LSCRIPT_CHECK_NAMED = $(call LSCRIPT,$2)
LSCRIPT_ENV_CHECK_NAMED = $(call LSCRIPT_ENV,$1,$3)
LSCRIPT_CHECK_NAMED_PARALLEL_MEM = $(call LSCRIPT_PARALLEL_MEM,$2,$3,$4,$5)
LSCRIPT_CHECK_NAMED_PARALLEL_LONG = $(call LSCRIPT_PARALLEL_LONG,$2,$3,$4,$5)
LSCRIPT_CHECK_NAMED_PARALLEL_SHORT = $(call LSCRIPT_PARALLEL_SHORT,$2,$3,$4,$5)
LSCRIPT_ENV_CHECK_NAMED_PARALLEL_MEM = $(call LSCRIPT_ENV_PARALLEL_MEM,$1,$3,$4,$5,$6)
LSCRIPT_ENV_CHECK_NAMED_PARALLEL_LONG = $(call LSCRIPT_ENV_PARALLEL_LONG,$1,$3,$4,$5,$6)
LSCRIPT_ENV_CHECK_NAMED_PARALLEL_SHORT = $(call LSCRIPT_ENV_PARALLEL_SHORT,$1,$3,$4,$5,$6)
LSCRIPT_CHECK_NAMED_MEM = $(call LSCRIPT_MEM,$2,$3,$4)
LSCRIPT_CHECK_NAMED_SHORT = $(call LSCRIPT_SHORT,$2,$3,$4)
LSCRIPT_CHECK_NAMED_LONG = $(call LSCRIPT_LONG,$2,$3,$4)
LSCRIPT_ENV_CHECK_NAMED_MEM = $(call LSCRIPT_ENV_MEM,$1,$3,$4,$5)
LSCRIPT_ENV_CHECK_NAMED_SHORT = $(call LSCRIPT_ENV_MEM,$1,$3,$4,$5)
LSCRIPT_ENV_CHECK_NAMED_LONG = $(call LSCRIPT_ENV_MEM,$1,$3,$4,$5)
else
#$(call LS,qsub,name,qsub_args,script)
LS = $(RUN_FIRST); (echo "umask $(UMASK); set -o pipefail; $(call ECHO_ACTIVATE_ENV,$(JRFLAB_MODULES_ENV));" $4) | $1 -N $(PROJECT_NAME)_$2 $3
#$(call LS,env,qsub,name,qsub_args,script)
LS_ENV = $(RUN_FIRST); (echo "umask $(UMASK); set -o pipefail; $(call ECHO_ACTIVATE_ENV,$1);" $5) | $2 -N $(PROJECT_NAME)_$3 $4
$(foreach wt,$(WALLTIMES),\
	$(eval LSCRIPT_$(wt) = $$(call LSCRIPT_NAMED_$(wt),$$(@F),$$1,$$2,$$3)) \
	$(eval LSCRIPT_ENV_$(wt) = $$(call LSCRIPT_ENV_NAMED_$(wt),$$1,$$(@F),$$2,$$3,$$4)) \
	$(eval LSCRIPT_PARALLEL_$(wt) = $$(call LSCRIPT_NAMED_PARALLEL_$(wt),$$(@F),$$1,$$2,$$3,$$4)) \
	$(eval LSCRIPT_ENV_PARALLEL_$(wt) = $$(call LSCRIPT_ENV_NAMED_PARALLEL_$(wt),$$1,$$(@F),$$2,$$3,$$4,$$5)) \
	$(eval LSCRIPT_NAMED_PARALLEL_$(wt) = $$(call LS,$$(QSUB),$$1,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_PARALLEL_$(wt),$$2,$$3,$$4),$$5)) \
	$(eval LSCRIPT_ENV_NAMED_PARALLEL_$(wt) = $$(call LS_ENV,$$1,$$(QSUB),$$2,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_PARALLEL_$(wt),$$3,$$4,$$5),$$6)) \
	$(eval LSCRIPT_NAMED_$(wt) = $$(call LS,$$(QSUB),$$1,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_$(wt),$$2,$$3),$$4)) \
	$(eval LSCRIPT_ENV_NAMED_$(wt) = $$(call LS_ENV,$$1,$$(QSUB),$$2,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_$(wt),$$3,$$4),$$5)) \
	$(eval LSCRIPT_CHECK_PARALLEL_$(wt) = $$(call LSCRIPT_CHECK_NAMED_PARALLEL_$(wt),$$(@F),$$1,$$2,$$3,$$4)) \
	$(eval LSCRIPT_ENV_CHECK_PARALLEL_$(wt) = $$(call LSCRIPT_ENV_CHECK_NAMED_PARALLEL_$(wt),$$1,$$(@F),$$2,$$3,$$4,$$5)) \
	$(eval LSCRIPT_CHECK_$(wt) = $$(call LSCRIPT_CHECK_NAMED_$(wt),$$(@F),$$1,$$2,$$3)) \
	$(eval LSCRIPT_ENV_CHECK_$(wt) = $$(call LSCRIPT_ENV_CHECK_NAMED_$(wt),$$1,$$(@F),$$2,$$3,$$4)) \
	$(eval LSCRIPT_DOCKER_PARALLEL_$(wt) = $$(call LSCRIPT_DOCKER_NAMED_PARALLEL_$(wt),$$(@F),$$1,$$2,$$3,$$4)) \
	$(eval LSCRIPT_DOCKER_$(wt) = $$(call LSCRIPT_DOCKER_NAMED_$(wt),$$(@F),$$1,$$2,$$3)) \
	$(eval LSCRIPT_CHECK_NAMED_PARALLEL_$(wt) = $$(call LS,$$(QSUB_CHECK),$$1,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_PARALLEL_$(wt),$$2,$$3,$$4),$$5)) \
	$(eval LSCRIPT_ENV_CHECK_NAMED_PARALLEL_$(wt) = $$(call LS_ENV,$$1,$$(QSUB_CHECK),$$3,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_PARALLEL_$(wt),$$3,$$4,$$5),$$6)) \
	$(eval LSCRIPT_DOCKER_NAMED_PARALLEL_$(wt) = $$(call LS,$$(QSUB),$$1,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_DOCKER_PARALLEL_$(wt),$$2,$$3,$$4),$$5)) \
	$(eval LSCRIPT_CHECK_NAMED_$(wt) = $$(call LS,$$(QSUB_CHECK),$$1,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_$(wt),$$2,$$3),$$4)) \
	$(eval LSCRIPT_ENV_CHECK_NAMED_$(wt) = $$(call LS_ENV,$$1,$$(QSUB_CHECK),$$2,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_$(wt),$$3,$$4),$$5)) \
	$(eval LSCRIPT_DOCKER_NAMED_$(wt) = $$(call LS,$$(QSUB_CHECK),$$1,$$($$(CLUSTER_ENGINE)_QSUB_ARGS) $$(call $$(CLUSTER_ENGINE)_DOCKER_$(wt),$$2,$$3),$$4)))

LSCRIPT = $(call LSCRIPT_NAMED_MEM,$(@F),1G,2G,$1)
LSCRIPT_ENV = $(call LSCRIPT_ENV_NAMED_MEM,$1,$(@F),1G,2G,$2)
LSCRIPT_NAMED = $(call LS,$(QSUB),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,1G,2G),$2)
LSCRIPT_ENV_NAMED = $(call LS_ENV,$1,$(QSUB),$2,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,1G,2G),$3)

LSCRIPT_CHECK = $(call LSCRIPT_CHECK_NAMED,$(@F),$1)
LSCRIPT_ENV_CHECK = $(call LSCRIPT_ENV_CHECK_NAMED,$1,$(@F),$2)

LSCRIPT_CHECK_NAMED = $(call LS,$(QSUB_CHECK),$1,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,1G,2G),$2)
LSCRIPT_ENV_CHECK_NAMED = $(call LS_ENV,$1,$(QSUB_CHECK),$2,$($(CLUSTER_ENGINE)_QSUB_ARGS) $(call $(CLUSTER_ENGINE)_MEM,1G,2G),$3)
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
CHECK_VCF = mkdir -p $(@D); if [ `grep -v '^\#' $< | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $< | wc -l` -eq 1 ]; then cp $< $@; else $1; fi
CHECK_MAF = mkdir -p $(@D); if [ `grep -v '^\#' $< | wc -l` -eq 1 ] && [ `grep '^Hugo_Symbol' $< | wc -l` -eq 1 ]; then cp $< $@; else $1; fi

# verify that vcf files have #CHROM
# usage: $(call VERIFY_VCF,tmp,final_dest)
VERIFY_VCF = if grep -q '^\#CHROM' $1; then mv $1 $2; else false; fi

<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
4< = $(word 4,$^)
<<<<< = $(word 5,$^)
5< = $(word 5,$^)

KNIT = $(RSCRIPT) modules/scripts/knit.R

PASS_FILTER_VCF = python modules/vcf_tools/pass_filter_vcf.py

define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef
ifeq ($(USE_CLUSTER),true)
ifeq ($(CLUSTER_ENGINE),SGE)
export DRMAA_LIBRARY_PATH = /common/sge/lib/lx24-amd64/libdrmaa.so.1.0
else
export DRMAA_LIBRARY_PATH = /opt/torque/lib/libdrmaa.so.1
endif
endif


endif
MAKEFILE_INC = true
