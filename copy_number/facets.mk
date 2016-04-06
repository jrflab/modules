# MSKCC facets R program
#
include modules/Makefile.inc

LOGDIR ?= log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

RUN_FACETS = $(RSCRIPT) modules/copy_number/runFacets.R
FACETS_PRE_CVAL ?= 50
FACETS_CVAL1 ?= 150
FACETS_CVAL2 ?= 50
FACETS_MIN_NHET ?= 25
FACETS_GATK_VARIANTS ?= false
FACETS_OPTS = --cval2 $(FACETS_CVAL2) --cval1 $(FACETS_CVAL1) --genome $(REF) --min_nhet $(FACETS_MIN_NHET) --pre_cval $(FACETS_PRE_CVAL)

GET_BASE_COUNTS = /ifs/e63data/reis-filho/usr/bin/GetBaseCounts
GET_BASE_COUNTS_OPTS = --filter_improper_pair --sort_output --maq 15 --baq 20 --cov 0 --fasta $(REF_FASTA)


ifeq ($(FACETS_GATK_VARIANTS),true)
FACETS_SNP_VCF = $(if $(TARGETS_FILE),facets/gatk_variant_input/all.variants.dbsnp.vcf.gz)
else
FACETS_SNP_VCF = $(if $(TARGETS_FILE),facets/targets_dbsnp.vcf.gz,$(DBSNP))
endif

MERGE_TN = $(PYTHON) /ifs/e63data/reis-filho/usr/bin/FACETS.app/mergeTN.py

FACETS_GENE_CN = $(RSCRIPT) modules/copy_number/facetsGeneCN.R
FACETS_GENE_CN_OPTS = $(if $(GENES_FILE),--genesFile $(GENES_FILE)) \
					  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) --mysqlDb $(EMBL_MYSQLDB_DB)

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/$(pair).cncf.txt) facets/geneCN.txt facets/geneCN_fill.txt

# FACETS_GATK_VARIANTS taget definitions
facets/gatk_variant_input/all.variants.dbsnp.vcf.gz : facets/gatk_variant_input/all.variants.snps.filtered.recode.vcf.gz
	$(INIT) $(BEDTOOLS) intersect -u -a $< -b $(DBSNP) | gzip -c > $@

facets/gatk_variant_input/all.variants.snps.filtered.recode.vcf.gz : $(foreach sample,$(UNMATCHED_SAMPLES),gatk/vcf/$(sample).variants.snps.filtered.vcf)
	$(RSCRIPT) modules/scripts/facets_gatk_variants.R

# no flag target definitions
facets/targets_dbsnp.vcf.gz : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -u -a $(DBSNP) -b $< | gzip -c > $@

facets/base_count/%.bc.gz : bam/%.bam $(FACETS_SNP_VCF) bam/%.bam.bai
	$(call LSCRIPT_CHECK_MEM,8G,13G,"$(GET_BASE_COUNTS) $(GET_BASE_COUNTS_OPTS) --bam $< --vcf $(<<) --out >( gzip -c > $@)")

define base-count-tumor-normal
facets/base_count/$1_$2.bc.gz : facets/base_count/$1.bc.gz facets/base_count/$2.bc.gz
	$$(call LSCRIPT_CHECK_MEM,8G,30G,"$$(MERGE_TN) $$^ | gzip -c > $$@")
endef

$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call base-count-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

facets/%.cncf.txt : facets/base_count/%.bc.gz
	$(call LSCRIPT_CHECK_MEM,8G,30G,"$(RUN_FACETS) $(FACETS_OPTS) --outPrefix facets/$* $<")

facets/geneCN.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/$(pair).cncf.txt)
	$(call LSCRIPT_CHECK_MEM,8G,30G,"$(FACETS_GENE_CN) $(FACETS_GENE_CN_OPTS) --outFile $@ $^")

facets/geneCN_fill.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/$(pair).cncf.txt) facets/geneCN.txt
	$(call LSCRIPT_CHECK_MEM,8G,30G,"$(FACETS_GENE_CN) $(FACETS_GENE_CN_OPTS) --outFile $@ $^")

include modules/bam_tools/processBam.mk
