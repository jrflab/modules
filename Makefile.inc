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
  SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_FILE))
endif

get_tumors = $(wordlist 1,$(shell expr $(words $1) - 1),$1)
get_normal = $(lastword $1)
get_space_sets = $(shell sed -n '/^\#/d; s/\s\+/ /g; s/\s\+$$//; $(1)p' $(SAMPLE_SET_FILE))
get_underscore_sets = $(shell sed -n '/^\#/d; s/\s\+/_/g; s/\s\+$$//; $(1)p' $(SAMPLE_SET_FILE))

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
  ifndef SAMPLE_SETS
    $(foreach i,$(SETS_SEQ), \
		$(eval SAMPLE_SETS += $(call get_underscore_sets,$i)))
  endif
  # tumor.SET => tumors
  $(foreach i,$(SETS_SEQ), \
	$(eval tumor.$(call get_underscore_sets,$i) := $(call get_tumors,$(set.$i))))
  # normal.SET => normal
  $(foreach i,$(SETS_SEQ), \
	$(eval normal.$(call get_underscore_sets,$i) := $(call get_normal,$(set.$i))))
  NORMAL_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_normal,$(set.$i)))
  TUMOR_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_tumors,$(set.$i)))
  SAMPLE_PAIRS ?= $(foreach set,$(SAMPLE_SETS),$(foreach tumor,$(tumor.$(set)),$(tumor)_$(normal.$(set))))
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

#ifeq ($(NORMAL_MATCHED),true)
#ifndef SAMPLE_PAIRS
#TUMOR_SAMPLES ?= $(shell sed 's/\w*$$//; s/ /\n/' $(SAMPLE_SET_FILE) | sort | uniq)
#NORMAL_SAMPLES ?= $(shell sed 's/.* //' $(SAMPLE_SET_FILE) | sort | uniq)
#$(foreach i,$(shell seq 1 $(NUM_SETS)),$(foreach tumor,$(wordlist 1,$(shell expr $(words $(set.$i)) - 1),$(set.$i)),$(eval normal_lookup.$(tumor) = $(lastword $(set.$i)))))
#$(foreach i,$(shell seq 1 $(NUM_SETS)),$(eval tumor_lookup.$(lastword $(set.$i)) += $(wordlist 1,$(shell expr $(words $(set.$i)) - 1),$(set.$i))))
#$(foreach tumor,$(TUMOR_SAMPLES),$(foreach normal,$(normal_lookup.$(tumor)),$(eval SAMPLE_PAIRS += $(tumor)_$(normal))))
#endif





#$(foreach i,$(shell seq 1 $(NUM_SETS)),$(info $(wordlist 1,$(shell expr $(words $(set.$i)) - 1),$(set.$i))))


ifneq ($(wildcard $(SAMPLE_SPLIT_FILE)),)
  A = $(shell cut -f1 $(SAMPLE_SPLIT_FILE))
  B = $(shell cut -f2 $(SAMPLE_SPLIT_FILE))
  $(foreach i,$(shell seq 1 $(words $(A))),$(eval split_lookup.$(word $i,$(A)) += $(word $i,$(B))))
  UNSPLIT_SAMPLES = $(B)
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
PATH := $(PATH):$(HOME)/share/usr/ensembl-git-tools/bin
export PATH
PYTHONPATH := $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7:/opt/common/python/python-2.7.5/lib:/opt/common/python/python-2.7.5/lib/python2.7:$(HOME)/share/usr/anaconda/lib/python2.7/site-packages:$(HOME)/share/usr/anaconda/lib/python2.7:$(HOME)/share/usr/anaconda/lib/python
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

PERL5LIB := /opt/common/Perl/lib/perl5/:$(HOME)/share/usr/perl5/lib/perl5:$(HOME)/share/usr/ActivePerl-5.18/lib
PERL5LIB := $(HOME)/share/usr/src/bioperl-live:$(PERL5LIB)
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-compara/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-variation/modules
PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-funcgen/modules

export PERL5LIB

#Interpreters
PYTHON = /opt/common/python/python-2.7.5/bin/python
ANACONDA_PYTHON = $(HOME)/share/usr/anaconda/bin/python
PERL = $(HOME)/share/usr/ActivePerl-5.18/bin/perl
#RSCRIPT = /usr/bin/Rscript
#R = /opt/common/R/R-3.0.0/bin/R
R = $(HOME)/share/usr/bin/R
#RSCRIPT = /opt/common/R/R-3.0.0/bin/Rscript
MY_RSCRIPT = $(HOME)/share/usr/bin/Rscript
RSCRIPT = $(MY_RSCRIPT)
#RSCRIPT = $(HOME)/usr/bin/Rscript

JARDIR := $(HOME)/share/usr/lib/java

### Applications
#BOWTIE = /opt/common/bowtie2/bowtie2-2.1.0/bowtie2
BOWTIE = $(HOME)/usr/bowtie2-2.2.4/bowtie2
#BWA = /opt/common/bwa/bwa-0.6.2-r126/bwa
BWA = $(HOME)/share/usr/bin/bwa-0.7.10
#PBWA = mpiexec $(HOME)/share/usr/bin/pBWA
#GSNAP = $(HOME)/share/usr/bin/gsnap
#GSNAP_REF = $(HOME)/share/references/gmapdb/
#GSNAP_MAPS = $(GSNAP_REF)/$(REF)/$(REF).maps
SAMTOOLS = $(HOME)/share/usr/bin/samtools-0.1.19
SAMTOOLS2 = $(HOME)/share/usr/bin/samtools
VCFUTILS = $(HOME)/share/usr/bin/vcfutils.pl
BCFTOOLS = /opt/bin/bcftools
BCFTOOLS2 = $(HOME)/share/usr/bin/bcftools
IGVTOOLS = $(HOME)/share/usr/IGVTools/igvtools
VCFTOOLS = $(HOME)/share/usr/bin/vcftools
VCF_SORT = $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl
SNP_EFF = $(JAVA) -Xmx4G -jar $(JARDIR)/snpEff.jar
CUFFLINKS = $(HOME)/share/usr/bin/cufflinks
CUFFCMP = $(HOME)/share/usr/bin/cuffcompare
#SNVMIX = $(HOME)/share/usr/bin/SNVMix2
#BAM2FASTQ = $(HOME)/share/usr/bin/bam2fastq
#QUALIMAP = $(HOME)/share/usr/qualimap/qualimap
#TOPHAT = $(HOME)/share/usr/tophat-2.0.8.Linux_x86_64/tophat2
#DEFUSE = /opt/common/defuse/defuse-0.6.1/scripts/defuse.pl
TMPDIR = $(HOME)/share/tmp

#VARSCAN
VARSCAN_JAR = $(JARDIR)/VarScan.v2.3.7.jar
VARSCAN_MEM = $(JAVA) -Xmx$1 -jar $(VARSCAN_JAR)
VARSCAN = $(call VARSCAN_MEM,8G)

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



# $(call strip-suffix,filename)
strip-suffix = $(firstword $(subst ., ,$(1)))


LINK = ln -svf $(notdir $1) $(2) && touch $1

NOW := $(shell date +"%F")


MKDIR = mkdir -p -m 775
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
INIT = $(MKDIRS); umask $(UMASK); $(LOG_MEM) ; set -o pipefail
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
QSUB_PRIORITY ?= -500
override QSUB_ARGS = -V -wd $(PWD) -now n -notify -o $(LOG) -j y -S $(SHELL) -p $(QSUB_PRIORITY) -b n
QMAKE_BIN = /common/sge/bin/lx24-amd64/qmake
QMAKE_QUEUE = $(QMAKE_BIN) -q $1
QMAKE = $(call QSUB_QUEUE,$(QUEUE))
QSUB_BIN = $(PERL) $(HOME)/share/scripts/qsub.pl $1 --
QSUB_QUEUE = $(call QSUB_BIN,$2) -q $1
QSUB = $(call QSUB_QUEUE,$(QUEUE))
# check output file
QSUB_CHECK = $(call QSUB_QUEUE,$(QUEUE),-o $@)
QSUB_CHECK_QUEUE = $(call QSUB_QUEUE,$1,-o $@)

USE_CLUSTER ?= true
ifeq ($(USE_CLUSTER),false)
LSCRIPT = $(INIT) echo "umask $(UMASK);" $1 | $(SHELL) > $(LOG).o 2> $(LOG).e
LSCRIPT_PARALLEL_MEM = $(INIT) echo "umask $(UMASK);" $4 | $(SHELL) > $(LOG).o 2> $(LOG).e
LSCRIPT_MEM = $(INIT) echo "umask $(UMASK);" $3 | $(SHELL) > $(LOG).o 2> $(LOG).e

LSCRIPT_NET = $(LSCRIPT)
LSCRIPT_PARALLEL_MEM_NET = $(LSCRIPT_PARALLEL_MEM)
LSCRIPT_MEM_NET = $(LSCRIPT_MEM)

LSCRIPT_NAMED = $(LSCRIPT)
LSCRIPT_NAMED_PARALLEL_MEM = $(LSCRIPT_PARALLEL_MEM)
LSCRIPT_NAMED_MEM = $(LSCRIPT_MEM)
else
LSCRIPT = $(call LSCRIPT_NAMED,$(@F),$1)
LSCRIPT_PARALLEL_MEM = $(call LSCRIPT_NAMED_PARALLEL_MEM,$(@F),$1,$2,$3,$4)
LSCRIPT_MEM = $(call LSCRIPT_NAMED_MEM,$(@F),$1,$2,$3)

LSCRIPT_NET = $(INIT) echo "umask $(UMASK);" $1 | $(call QSUB_QUEUE,out.q) -N X$(@F) $(QSUB_ARGS) -l internet=1
LSCRIPT_BIG = $(INIT) echo "umask $(UMASK);" $1 | $(call QSUB_QUEUE,big.q) -N X$(@F) $(QSUB_ARGS)
LSCRIPT_PARALLEL_MEM_NET = $(INIT) echo "umask $(UMASK);" $4 | $(call QSUB_QUEUE,out.q) -N X$(@F) $(QSUB_ARGS) -pe $(PARALLEL_ENV) $1 $(call MEM_FREE,$2,$3) -l internet=1
LSCRIPT_PARALLEL_MEM_BIG = $(INIT) echo "umask $(UMASK);" $4 | $(call QSUB_QUEUE,big.q) -N X$(@F) $(QSUB_ARGS) -pe $(PARALLEL_ENV) $1 $(call MEM_FREE,$2,$3)
LSCRIPT_MEM_NET = $(INIT) echo "umask $(UMASK);" $3 | $(call QSUB_QUEUE,out.q) -N X$(@F) $(QSUB_ARGS) $(call MEM_FREE,$1,$2) -l internet=1
LSCRIPT_MEM_BIG = $(INIT) echo "umask $(UMASK);" $3 | $(call QSUB_QUEUE,big.q) -N X$(@F) $(QSUB_ARGS) $(call MEM_FREE,$1,$2)

LSCRIPT_NAMED = $(INIT) echo "umask $(UMASK);" $2 | $(QSUB) -N X$1 $(QSUB_ARGS)
LSCRIPT_NAMED_PARALLEL_MEM = $(INIT) echo "umask $(UMASK);" $5 | $(QSUB) -N X$1 $(QSUB_ARGS) -pe $(PARALLEL_ENV) $2 $(call MEM_FREE,$3,$4)
LSCRIPT_NAMED_MEM = $(INIT) echo "umask $(UMASK);" $4 | $(QSUB) -N X$1 $(QSUB_ARGS) $(call MEM_FREE,$2,$3)

LSCRIPT_CHECK = $(call LSCRIPT_CHECK_NAMED,$(@F),$1)
LSCRIPT_CHECK_PARALLEL_MEM = $(call LSCRIPT_CHECK_NAMED_PARALLEL_MEM,$(@F),$1,$2,$3,$4)
LSCRIPT_CHECK_MEM = $(call LSCRIPT_CHECK_NAMED_MEM,$(@F),$1,$2,$3)

LSCRIPT_CHECK_NET = $(INIT) echo "umask $(UMASK);" $1 | $(call QSUB_CHECK_QUEUE,out.q) -N X$(@F) $(QSUB_ARGS) -l internet=1
LSCRIPT_CHECK_PARALLEL_MEM_NET = $(INIT) echo "umask $(UMASK);" $4 | $(call QSUB_CHECK_QUEUE,out.q) -N X$(@F) $(QSUB_ARGS) -pe $(PARALLEL_ENV) $1 $(call MEM_FREE,$2,$3) -l internet=1
LSCRIPT_CHECK_MEM_NET = $(INIT) echo "umask $(UMASK);" $3 | $(call QSUB_CHECK_QUEUE,out.q) -N X$(@F) $(QSUB_ARGS) $(call MEM_FREE,$1,$2) -l internet=1

LSCRIPT_CHECK_NAMED = $(INIT) echo "umask $(UMASK);" $2 | $(QSUB_CHECK) -N X$1 $(QSUB_ARGS)
LSCRIPT_CHECK_NAMED_PARALLEL_MEM = $(INIT) echo "umask $(UMASK);" $5 | $(QSUB_CHECK) -N X$1 $(QSUB_ARGS) -pe $(PARALLEL_ENV) $2 $(call MEM_FREE,$3,$4)
LSCRIPT_CHECK_NAMED_MEM = $(INIT) echo "umask $(UMASK);" $4 | $(QSUB_CHECK) -N X$1 $(QSUB_ARGS) $(call MEM_FREE,$2,$3)
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

SOMATIC_FILTER_VCF = $(PERL) $(HOME)/share/scripts/somaticFilterVCF.pl

GET_INSERT_SIZE = $(HOME)/share/usr/bin/getInsertSize.py

<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
4< = $(word 4,$^)
<<<<< = $(word 5,$^)
5< = $(word 5,$^)

#signature md5
#makepp_no_builtin = 1
#makepp_percent_subdirs := 1

define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef
