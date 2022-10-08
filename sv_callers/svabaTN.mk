include modules/Makefile.inc

LOGDIR = log/svabaTN.$(NOW)

SVABA_CORES ?= 8
SVABA_MEM_CORE ?= 6G
SVABA_REF ?= $(REF_FASTA)
SVABA_DBSNP ?= $(HOME)/share/lib/resource_files/svaba/dbsnp_indel.vcf
SVABA_BLACKLIST ?= $(HOME)/share/lib/resource_files/wgs_blacklist_meres.bed
SVABA_ENV ?= $(HOME)/share/usr/env/svaba-1.1.0
SVABA ?= svaba

svaba : $(foreach pair,$(SAMPLE_PAIRS),svaba/$(pair).svaba.somatic.indel.vcf)

define svaba-tumor-normal
svaba/$1_$2.svaba.somatic.indel.vcf : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -n $(SVABA_CORES) -s 4G -m $(SVABA_MEM_CORE) -v $(SVABA_ENV) -w 72:00:00,"set -o pipefail && \
												 mkdir -p svaba && \
										 		 cd svaba && \
												 $$(SVABA) run \
												 -t ../bam/$1.bam \
												 -n ../bam/$2.bam \
												 -p $$(SVABA_CORES) \
												 -D $$(SVABA_DBSNP) \
												 -L 100000 \
												 -x 25000 \
												 -k $$(SVABA_BLACKLIST) \
												 -a $1_$2 \
												 -G $$(SVABA_REF)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call svaba-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


..DUMMY := $(shell mkdir -p version; \
	     $(SVABA) --help &> version/svabaTN.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: svaba
