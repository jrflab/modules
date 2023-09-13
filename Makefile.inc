ifneq ("$(wildcard project_config.inc)", "")
	include project_config.inc
endif
include modules/config.inc
ifdef PROJECT_CONFIG
include $(PROJECT_CONFIG)
endif
ifneq ("$(wildcard config.inc)", "")
	include config.inc
endif

REF ?= b37
include modules/genome_inc/$(REF).inc

ifndef MAKEFILE_INC

SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt

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
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set_lookup.$(sample) := $i)))
  $(foreach i,$(SETS_SEQ), \
	  $(eval SAMPLE_SETS += $(call get_underscore_sets,$i)))
  $(foreach i,$(SETS_SEQ), \
	$(eval tumor.$(call get_underscore_sets,$i) := $(call get_tumors,$(set.$i))))
  $(foreach i,$(SETS_SEQ), \
	$(eval normal.$(call get_underscore_sets,$i) := $(call get_normal,$(set.$i))))
  NORMAL_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_normal,$(set.$i)))
  TUMOR_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_tumors,$(set.$i)))
  SAMPLE_PAIRS = $(foreach set,$(SAMPLE_SETS),$(foreach tumor,$(tumor.$(set)),$(tumor)_$(normal.$(set))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor.$(tumor)_$(normal.$(set)) := $(tumor))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval normal.$(tumor)_$(normal.$(set)) := $(normal.$(set)))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(tumor) := $(tumor)_$(normal.$(set)))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(normal.$(set)) := $(tumor)_$(normal.$(set)))))
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

export TMPDIR := $(HOME)/share/data/$(USER)/tmp

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
COLLECT_ALIGNMENT_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectAlignmentSummaryMetrics $(PICARD_OPTS)
COLLECT_INSERT_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectInsertSizeMetrics $(PICARD_OPTS)
COLLECT_OXOG_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectOxoGMetrics $(PICAD_OPTS)
COLLECT_GC_BIAS = $(PICARD) -Xmx$(PICARD_MEM) CollectGcBiasMetrics $(PICARD_OPTS)
BAM_INDEX = $(PICARD) -Xmx$(PICARD_MEM) BamIndexStats $(PICARD_OPTS)
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

SNP_EFF = $(call SNP_EFF_MEM,$(DEFAULT_JAVA_MEM))
SNP_EFF_MEM = $(JAVA8) -Xmx$1 -jar $(SNP_EFF_JAR)
SNP_SIFT_MEM = $(JAVA8) -Xmx$1 -jar $(SNP_SIFT_JAR)
SNP_SIFT = $(call SNP_SIFT_MEM,$(DEFAULT_JAVA_MEM))
VCF_EFF_ONE_PER_LINE = $(HOME)/share/usr/snpEff-4.1/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = modules/scripts/join_eff.pl

COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)

JAVA_ARGS = -Djava.io.tmpdir=$(TMPDIR)
JAVA := java $(JAVA_ARGS)
JAVA6 := $(JAVA6_BIN) $(JAVA_ARGS)
JAVA7 := $(JAVA7_BIN) $(JAVA_ARGS)
JAVA8 := $(JAVA8_BIN) $(JAVA_ARGS)

NO_RM ?= false
ifeq ($(NO_RM),true)
    RM := touch
    RMR = touch
else
    RM := rm -f
    RMR := rm -r
endif

, := ,
space :=
space +=
$(space) := 
$(space) +=

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

### SGE
QUEUE ?=
WALLTIMES = MEM SHORT LONG
SGE_PARALLEL_ENV = smp
SGE_NODES = e02 e06
SGE_QSUB_SCRIPT = python modules/scripts/qsub.py
export SGE_ROOT = /common/sge

### PBS
PBS_NODES = gpu-2-14 cpu-6-1 gpu-1-4

MAKE_PID := $(shell echo $$PPID)
JOB_FLAG := $(filter -j%, $(subst -j ,-j,$(shell ps T | grep "^\s*$(MAKE_PID).*$(MAKE)")))
JOBS := $(subst -j,,$(JOB_FLAG))

### LOCAL
RUN_FIRST = umask $(UMASK); $(if $(and $(findstring true,$(USE_CLUSTER)),$(JOBS)),sleep $$[ ( $$RANDOM % 20 ) + 1 ]s;) $(MKDIRS)
INIT = $(RUN_FIRST); $(call ACTIVATE_ENV,$(JRFLAB_MODULES_ENV)); umask $(UMASK); set -o pipefail;
INIT_ENV = $(RUN_FIRST); $(call ACTIVATE_ENV,$1); umask $(UMASK); set -o pipefail;

CLUSTER_ENGINE ?= SGE
USE_CLUSTER ?= true

RUN_OPTS = $(if $(findstring false,$(USE_CLUSTER)),-l) -o $@ $(if $(CLUSTER_ENGINE),-g $(CLUSTER_ENGINE)) \
		   --default_env $(JRFLAB_MODULES_ENV) -S $(SHELL) \
		   --log_file $(LOG) $(if $($(CLUSTER_ENGINE)_NODES),--servers $($(CLUSTER_ENGINE)_NODES)) \
		   $(if $(PROJECT_NAME),-p $(PROJECT_NAME))
RUN = echo $2 | python modules/scripts/run.py $(RUN_OPTS) $1

MD5 = md5sum $(@:.md5=) > $@
CHECK_MD5 = for i in {0..20}; do if md5sum -c $(filter %.md5,$^); then break; fi; sleep 5; done;
<M = $(<:.md5=)
^M = $(^:.md5=)
@M = $(@:.md5=)
<<M = $(word 2,$(^M))
<<<M = $(word 3,$(^M))

# check for valid header i.e. non-empty vcf file
# usage: $(call CHECK_VCF,vcf.file,target,command)
CHECK_VCF = mkdir -p $(@D); if [ `grep -v '^\#' $< | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $< | wc -l` -eq 1 ]; then cp $< $@; else $1; fi
CHECK_MAF = mkdir -p $(@D); if [ `grep -v '^\#' $< | wc -l` -eq 1 ] && [ `grep '^Hugo_Symbol' $< | wc -l` -eq 1 ]; then cp $< $@; else $1; fi
CHECK_UVCF = mkdir -p $(@D); if [ `grep -v '^\#' $< | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $< | wc -l` -eq 1 ]; then sed 's/FILTER\tINFO/FILTER\tUPS-COORDINATE\tINFO/' $< > $@; else $1; fi

# verify that vcf files have #CHROM
# usage: $(call VERIFY_VCF,tmp,final_dest)
VERIFY_VCF = if grep -q '^\#CHROM' $1; then mv $1 $2 && if [ -s $1.idx ]; then mv $1.idx $2.idx; fi ; else false; fi
VERIFY_VCF_GZ = if zgrep -q '^\#CHROM' $1; then mv $1 $2 && if [ -s $1.tbi ]; then mv $1.tbi $2.tbi; fi ; else false; fi

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
ifeq ($(CLUSTER_ENGINE),LSF)
export DRMAA_LIBRARY_PATH = /admin/lsf/10.1/linux3.10-glibc2.17-x86_64/lib/libdrmaa.so
else
export DRMAA_LIBRARY_PATH = /opt/torque/lib/libdrmaa.so.1
endif
endif
endif

endif
MAKEFILE_INC = true
