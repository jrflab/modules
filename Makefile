ifneq ("$(wildcard config.inc)", "")
	include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
	include project_config.inc
endif
include modules/config/config.inc

export

NUM_ATTEMPTS ?= 10
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

USE_CLUSTER ?= true
QMAKE = modules/scripts/runtime/qmake.pl -n $@.$(NOW) $(if $(SLACK_CHANNEL),-c $(SLACK_CHANNEL)) -r $(NUM_ATTEMPTS) -m -s -- make
NUM_JOBS ?= 250

define RUN_QMAKE
$(QMAKE) -e -f $1 -j $2 $(TARGET) && \
	mkdir -p completed_tasks && \
	touch completed_tasks/$@
endef

RUN_MAKE = $(if $(findstring false,$(USE_CLUSTER))$(findstring n,$(MAKEFLAGS)),+$(MAKE) -f $1,$(call RUN_QMAKE,$1,$(NUM_JOBS)))

#==================================================
# FASTQ / BAM file aligners
#==================================================

TARGETS += bwa_mem
bwa_mem :
	$(call RUN_MAKE,modules/fastq_aligners/bwa_mem.mk)
	
TARGETS += bwa_wgs
bwa_wgs :
	$(call RUN_MAKE,modules/fastq_aligners/bwa_wgs.mk)
	
#==================================================
# BAM file utilities
#==================================================

TARGETS += split_rg
splite_rg :
	$(call RUN_MAKE,modules/bam_tools/split_rg.mk)

#==================================================
# FASTQ file utilities
#==================================================

TARGETS += subsample_fastq
subsample_fastq :
	$(call RUN_MAKE,modules/fastq_tools/subsample_fastq.mk)

#==================================================
# VCF file utilities
#==================================================

TARGETS += annotate_vcf_context
annotate_vcf_context :
	$(call RUN_MAKE,modules/vcf_tools/annotate_vcf_context.mk)
	
TARGETS += annotate_vcf_maf
annotate_vcf_maf :
	$(call RUN_MAKE,modules/vcf_tools/annotate_vcf_maf.mk)
	
TARGETS += annotate_maf_vcf
annotate_maf_vcf :
	$(call RUN_MAKE,modules/vcf_tools/annotate_maf_vcf.mk)

#==================================================
# BETA testing
#==================================================

TARGETS += pileup_metrics
pileup_metrics :
	$(call RUN_MAKE,modules/mission_bio/pileup_metrics.mk)
	
TARGETS += get_basecount
get_basecount :
	$(call RUN_MAKE,modules/mission_bio/get_basecount.mk)
	

.PHONY : $(TARGETS)
