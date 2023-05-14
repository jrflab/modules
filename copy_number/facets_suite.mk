include modules/Makefile.inc

LOGDIR ?= log/facets_suite.$(NOW)

FACETS_MAX_DEPTH ?= 15000
FACETS_PRE_CVAL ?= 50
FACETS_CVAL1 ?= 150
FACETS_CVAL2 ?= 50
FACETS_MIN_NHET ?= 25
FACETS_SNP_NBHD ?= 250
FACETS_HET_THRESHOLD ?= 0.25

facets_suite : facets_suite/vcf/targets_dbsnp.vcf
#	       $(foreach pair,$(SAMPLE_PAIRS),facets_suite/$(pair)/$(pair).snp_pileup.gz)

facets_suite/vcf/targets_dbsnp.vcf : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@
	

#define snp-pileup
#facets_suite/$1_$2/$1_$2.snp_pileup.gz : facets_suite/vcf/targets_dbsnp.vcf bam/$1.bam bam/$2.bam
#	$$(call RUN,-c -s 1G -m 2G -v $(FACETS_SUITE_ENV),"set -o pipefail && \
#							   snp-pileup-wrapper.R --verbose \
#							   -sp /home/$(USER)/share/usr/env/r-facets-suite-2.0.8/bin/snp-pileup \
#							   --vcf-file $$(<) \
#							   --tumor-bam $$(<<) \
#							   --normal-bam $$(<<<) \
#							   --output-prefix facets_suite/$1_$2/$1_$2 \
#							   --pseudo-snps NULL \
#							   --max-depth $$(FACETS_MAX_DEPTH)")
#	
#endef
#$(foreach pair,$(SAMPLE_PAIRS),\
#	$(eval $(call snp-pileup,\
#		$(tumor.$(pair)),$(normal.$(pair)))))



..DUMMY := $(shell mkdir -p version; \
	     $(FACETS_SUITE_ENV)/bin/R --version > version/facets_suite.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: facets_suite
