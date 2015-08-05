# MSKCC facets R program
#
include modules/Makefile.inc

LOGDIR = log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

RUN_FACETS = $(RSCRIPT) modules/copy_number/runFacets.R
FACETS_PRE_CVAL ?= 50
FACETS_CVAL1 ?= 150
FACETS_CVAL2 ?= 50
FACETS_MIN_NHET ?= 25
FACETS_OPTS = --cval2 $(FACETS_CVAL2) --cval1 $(FACETS_CVAL1) --genome $(REF) --min_nhet $(FACETS_MIN_NHET) --pre_cval $(FACETS_PRE_CVAL)

GET_BASE_COUNTS = /ifs/e63data/socci/Code/FACETS/bin/GetBaseCounts
GET_BASE_COUNTS_OPTS = --filter_improper_pair --sort_output --maq 15 --baq 20 --cov 0 --fasta $(REF_FASTA)

FACETS_SNP_VCF = $(if $(TARGETS_FILE),facets/targets_dbsnp.vcf.gz,$(DBSNP))

MERGE_TN = $(PYTHON) /ifs/e63data/socci/Code/FACETS/mergeTN.py

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/$(pair).out)

facets/targets_dbsnp.vcf.gz : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -u -a $(DBSNP) -b $< | gzip -c > $@

facets/base_count/%.bc.gz : bam/%.bam $(FACETS_SNP_VCF)
	$(call LSCRIPT_MEM,8G,10G,"$(GET_BASE_COUNTS) $(GET_BASE_COUNTS_OPTS) --bam $< --vcf $(<<) --out >( gzip -c > $@)")

define base-count-tumor-normal
facets/base_count/$1_$2.bc.gz : facets/base_count/$1.bc.gz facets/base_count/$2.bc.gz
	$$(call LSCRIPT_MEM,8G,10G,"$$(MERGE_TN) $$^ | gzip -c > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call base-count-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

facets/%.out : facets/base_count/%.bc.gz
	$(call LSCRIPT_MEM,8G,10G,"$(RUN_FACETS) $(FACETS_OPTS) --outPrefix facets/$* $<")

