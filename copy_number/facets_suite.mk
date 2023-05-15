include modules/Makefile.inc

LOGDIR ?= log/facets_suite.$(NOW)

FACETS_MAX_DEPTH ?= 15000
FACETS_CVAL ?= 50
FACETS_PURITY_CVAL ?= 30
FACETS_MIN_NHET ?= 15
FACETS_PURITY_MIN_NHET ?= 10
SNP_WINDOW_SIZE ?= 250
NORMAL_DEPTH ?= 25

facets_suite : facets_suite/vcf/targets_dbsnp.vcf \
	       $(foreach pair,$(SAMPLE_PAIRS),facets_suite/$(pair)/$(pair).snp_pileup.gz) \
	       $(foreach pair,$(SAMPLE_PAIRS),facets_suite/$(pair)/taskcomplete) \
	       facets_suite/summary/summary.txt

facets_suite/vcf/targets_dbsnp.vcf : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@
	

define snp-pileup
facets_suite/$1_$2/$1_$2.snp_pileup.gz : facets_suite/vcf/targets_dbsnp.vcf bam/$1.bam bam/$2.bam
	$$(call RUN,-c -s 2G -m 4G -v $(FACETS_SUITE_ENV),"set -o pipefail && \
							   snp-pileup-wrapper.R --verbose \
							   -sp /home/$(USER)/share/usr/env/r-facets-suite-2.0.8/bin/snp-pileup \
							   --vcf-file $$(<) \
							   --tumor-bam $$(<<) \
							   --normal-bam $$(<<<) \
							   --output-prefix facets_suite/$1_$2/$1_$2 \
							   --pseudo-snps NULL \
							   --max-depth $$(FACETS_MAX_DEPTH)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call snp-pileup,$(tumor.$(pair)),$(normal.$(pair)))))

define run-facets
facets_suite/$1_$2/taskcomplete : facets_suite/$1_$2/$1_$2.snp_pileup.gz
	$$(call RUN,-c -s 4G -m 6G -v $(FACETS_SUITE_ENV),"set -o pipefail && \
							   run-facets-wrapper.R --verbose \
							   --counts-file $$(<) \
							   --sample-id $1_$2 \
							   --directory facets_suite/$1_$2/ \
							   --everything \
							   --genome hg19 \
							   --cval $$(FACETS_CVAL) \
							   --purity-cval $$(FACETS_PURITY_CVAL) \
							   --min-nhet $$(FACETS_MIN_NHET) \
							   --purity-min-nhet $$(FACETS_PURITY_MIN_NHET) \
							   --snp-window-size $$(SNP_WINDOW_SIZE) \
							   --normal-depth $$(NORMAL_DEPTH) \
							   --seed 0 \
							   --legacy-output True \
							   --facets-lib-path /home/$(USER)/share/usr/env/r-facets-suite-2.0.8/lib/R/library/ && \
							   echo 'finished!' > $$(@)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call run-facets,$(tumor.$(pair)),$(normal.$(pair)))))
		

facets_suite/summary/summary.txt : $(foreach pair,$(SAMPLE_PAIRS),facets_suite/$(pair)/taskcomplete)
	$(call RUN, -c -n 1 -s 24G -m 48G -v $(INNOVATION_ENV),"set -o pipefail && \
								$(RSCRIPT) $(SCRIPTS_DIR)/facets_suite.R --option 1 --sample_pairs '$(SAMPLE_PAIRS)'")
					  

..DUMMY := $(shell mkdir -p version; \
	     $(FACETS_SUITE_ENV)/bin/R --version > version/facets_suite.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: facets_suite
