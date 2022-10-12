ifneq ("$(wildcard project_config.inc)", "")
	include project_config.inc
endif
include modules/config/config.inc
ifdef PROJECT_CONFIG
include $(PROJECT_CONFIG)
endif
ifneq ("$(wildcard config.inc)", "")
	include config.inc
endif

include modules/genome_inc/$(REF).inc

ifndef MAKEFILE_INC

SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt

PROJECT_NAME = $(shell pwd | sed 's:.*/home/::; s:.*/share/::; s:.*/$(USER)/::; s:.*/juno/::; s:.*/work/::; s:.*/reis-filho/::; s:/:_:g')

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

COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)

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
MKDIR = mkdir -p -m 775
MKDIRS = $(MKDIR) $(LOGDIR)/$(@D) $(@D)
LOG = $(PWD)/$(LOGDIR)/$(@).log
NOW := $(shell date +"%F")

UMASK = 002

export PATH := $(INNOVATION_LAB_ENV)/bin:$(PATH)
ECHO_ACTIVATE_ENV = unset PYTHONPATH; OLDPATH=\$$PATH; tries=0; until [[ \$$tries -gt 10 ]] || source $1/bin/activate $1; do export PATH=\$$OLDPATH; tries=\`expr \$$tries + 1\`; jitter=\`expr \$$RANDOM % 10 + 1\`; sleep \`expr \$$jitter \\* \$$tries\`; done
ACTIVATE_ENV = unset PYTHONPATH; OLDPATH=$$PATH; tries=0; until [[ $$tries -gt 10 ]] || source $1/bin/activate $1; do export PATH=$$OLDPATH; tries=`expr $$tries + 1`; jitter=`expr $$RANDOM % 10 + 1`; sleep `expr $$jitter \* $$tries`; done

MAKE_PID := $(shell echo $$PPID)
JOB_FLAG := $(filter -j%, $(subst -j ,-j,$(shell ps T | grep "^\s*$(MAKE_PID).*$(MAKE)")))
JOBS := $(subst -j,,$(JOB_FLAG))

RUN_FIRST = umask $(UMASK); $(if $(and $(findstring true,$(USE_CLUSTER)),$(JOBS)),sleep $$[ ( $$RANDOM % 20 ) + 1 ]s;) $(MKDIRS)
INIT = $(RUN_FIRST); $(call ACTIVATE_ENV,$(INNOVATION_LAB_ENV)); umask $(UMASK); set -o pipefail;
INIT_ENV = $(RUN_FIRST); $(call ACTIVATE_ENV,$1); umask $(UMASK); set -o pipefail;

USE_CLUSTER ?= true
CLUSTER_ENGINE ?= LSF

RUN_OPTS = $(if $(findstring false,$(USE_CLUSTER)),-l) -o $@ $(if $(CLUSTER_ENGINE),-g $(CLUSTER_ENGINE)) \
		   --default_env $(INNOVATION_LAB_ENV) -S $(SHELL) \
		   --log_file $(LOG) $(if $($(CLUSTER_ENGINE)_NODES),--servers $($(CLUSTER_ENGINE)_NODES)) \
		   $(if $(PROJECT_NAME),-p $(PROJECT_NAME))
RUN = echo $2 | python $(SCRIPTS_DIR)/runtime/run.py $(RUN_OPTS) $1

MD5 = md5sum $(@:.md5=) > $@
CHECK_MD5 = for i in {0..20}; do if md5sum -c $(filter %.md5,$^); then break; fi; sleep 5; done;
<M = $(<:.md5=)
^M = $(^:.md5=)
@M = $(@:.md5=)
<<M = $(word 2,$(^M))
<<<M = $(word 3,$(^M))
<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
4< = $(word 4,$^)
<<<<< = $(word 5,$^)
5< = $(word 5,$^)

define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef

ifeq ($(USE_CLUSTER),true)
ifeq ($(CLUSTER_ENGINE),SGE)
export DRMAA_LIBRARY_PATH = /common/sge/lib/lx24-amd64/libdrmaa.so
else
ifeq ($(CLUSTER_ENGINE),LSF)
export DRMAA_LIBRARY_PATH = /admin/lsfjuno/lsf/10.1/linux3.10-glibc2.17-x86_64/lib/libdrmaa.so
endif
endif
endif

endif
MAKEFILE_INC = true
