# MSKCC facets R program
#
include modules/Makefile.inc

LOGDIR ?= log/facets.$(NOW)

#.SECONDARY:
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
NORMAL_BASE_COUNT_CENTER = 100
NORMAL_BASE_COUNT_CLIP_THRESHOLD = 800
RECENTER_BASE_COUNTS = python modules/copy_number/recenter_base_count.py
RECENTER_BASE_COUNTS_OPTS = --normal_center $(NORMAL_BASE_COUNT_CENTER) --threshold $(NORMAL_BASE_COUNT_CLIP_THRESHOLD)

FACETS_DBSNP = $(if $(TARGETS_FILE),facets/vcf/targets_dbsnp.vcf.gz,$(DBSNP))

# downsample bams using GATK (max depth = 800)
FACETS_DOWNSAMPLE_BAM ?= false
# recenter and clip the normal base counts
FACETS_RECENTER_CLIP_BASE_COUNTS ?= false

# augment dbsnp with calls from heterozygous calls from gatk
FACETS_UNION_GATK_DBSNP ?= false
ifeq ($(FACETS_UNION_GATK_DBSNP),true)
FACETS_SNP_VCF = facets/vcf/dbsnp_het_gatk.snps.vcf.gz
else
FACETS_SNP_VCF = $(FACETS_DBSNP)
endif

MERGE_TN = $(PYTHON) /ifs/e63data/reis-filho/usr/bin/FACETS.app/mergeTN.py

FACETS_GENE_CN = $(RSCRIPT) modules/copy_number/facetsGeneCN.R
FACETS_FILL_GENE_CN = $(RSCRIPT) modules/copy_number/facetsFillGeneCN.R
FACETS_GENE_CN_OPTS = $(if $(GENES_FILE),--genesFile $(GENES_FILE)) \
					  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) \
					  --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) \
					  --mysqlDb $(EMBL_MYSQLDB_DB)
FACETS_PLOT_GENE_CN = $(RSCRIPT) modules/copy_number/facetsGeneCNPlot.R
FACETS_PLOT_GENE_CN_OPTS = --sampleColumnPostFix '_LRR_threshold'


#------
# rules
#------

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt) \
	facets/geneCN.txt facets/geneCN.pdf facets/geneCN.fill.txt

facets/vcf/dbsnp_het_gatk.snps.vcf : $(FACETS_DBSNP:.gz=) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.het.pass.vcf)
	$(call LSCRIPT_CHECK_MEM,4G,6G,"$(call GATK_MEM,3G) -T CombineVariants --minimalVCF $(foreach i,$^, --variant $i) -R $(REF_FASTA) -o $@")

# flag homozygous calls
%.het.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call GATK_MEM,8G) -V $< -T VariantFiltration -R $(REF_FASTA) --genotypeFilterName 'hom' --genotypeFilterExpression 'isHet == 0' -o $@")

%.vcf.gz : %.vcf
	$(INIT) cat $< | gzip -c > $@

%.vcf : %.vcf.gz
	$(INIT) zcat $< > $@

# no flag target definitions
facets/vcf/targets_dbsnp.vcf.gz : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< | gzip -c > $@

facets/base_count/%.bc.gz : $(if $(findstring $(FACETS_DOWNSAMPLE_BAM),true),bam/%.dcov800.bam bam/%.dcov800.bam.bai,bam/%.bam bam/%.bam.bai) $(FACETS_SNP_VCF)
	$(call LSCRIPT_CHECK_MEM,8G,13G,"$(GET_BASE_COUNTS) $(GET_BASE_COUNTS_OPTS) --bam $< --vcf $(<<<) --out >( gzip -c > $@)")

define base-count-tumor-normal
facets/base_count/$1_$2.bc.gz : facets/base_count/$1.bc.gz facets/base_count/$2.bc.gz
	$$(call LSCRIPT_CHECK_MEM,8G,30G,"$$(MERGE_TN) $$^ \
		$$(if $$(findstring true,$$(FACETS_RECENTER_CLIP_BASE_COUNTS)), \
		| $$(RECENTER_BASE_COUNTS) $$(RECENTER_BASE_COUNTS_OPTS)) \
		| gzip -c > $$@")
endef

$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call base-count-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

facets/cncf/%.cncf.txt : facets/base_count/%.bc.gz
	$(call LSCRIPT_CHECK_MEM,8G,30G,"$(RUN_FACETS) $(FACETS_OPTS) --outPrefix $(@D)/$* $<")

facets/geneCN.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call LSCRIPT_CHECK_MEM,8G,30G,"$(FACETS_GENE_CN) $(FACETS_GENE_CN_OPTS) --outFile $@ $^")


#-------------------
# facets fill script
#-------------------

GENECN_TXT = facets/geneCN.txt
GENECN_PDF = facets/geneCN.raw.pdf
GENECN_FILL_TXT = facets/geneCN.fill.txt
GENECN_FILL_PDF = facets/geneCN.fill.pdf

facets/geneCN.fill.txt : $(GENECN_TXT) $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call LSCRIPT_CHECK_MEM,8G,30G,"$(FACETS_FILL_GENE_CN) \
		--geneCN_fill_pdf $(GENECN_FILL_PDF) \
		--geneCN_fill_txt $(GENECN_FILL_TXT) \
		--geneCN_pdf $(GENECN_PDF) \
		--geneCN_txt $(GENECN_TXT) \
		$(filter %.cncf.txt,$^)")


#----------------
# geneCN plotting
#----------------

facets/geneCN.pdf : facets/geneCN.txt
	$(call LSCRIPT_MEM,8G,10G,"$(FACETS_PLOT_GENE_CN) $(FACETS_PLOT_GENE_CN_OPTS) $< $@")


#---------
# includes
#---------

include modules/variant_callers/gatk.mk
include modules/bam_tools/processBam.mk
