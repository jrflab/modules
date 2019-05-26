include modules/Makefile.inc

LOGDIR ?= log/fastaccess.$(NOW)
PHONY += fastaccess fastaccess/vcf fastaccess/pileup fastaccess/cncf fastaccess/plots

FACETS_ENV = $(HOME)/share/usr/anaconda-envs/facets-0.5.6/
RUN_FACETS = $(RSCRIPT) modules/copy_number/runFacets.R
PLOT_FACETS = $(RSCRIPT) modules/copy_number/plotFacets.R
POOL_A_BED = ~/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
POOL_B_BED = ~/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
FACETS_PRE_CVAL ?= 50
FACETS_CVAL1 ?= 150
FACETS_CVAL2 ?= 50
FACETS_MIN_NHET ?= 25
FACETS_SNP_NBHD ?= 250
FACETS_NDEPTH_MAX ?= 100000
FACETS_HET_THRESHOLD ?= 0.25
FACETS_OPTS = --genome $(REF) \
			  $(if $(facets_cval1.$1),--cval1 $(facets_diplogr.$1),--cval1 $(FACETS_CVAL1)) \
			  $(if $(facets_cval2.$1),--cval2 $(facets_diplogr.$1),--cval2 $(FACETS_CVAL1)) \
			  --het_threshold $(FACETS_HET_THRESHOLD) \
			  --min_nhet $(FACETS_MIN_NHET) \
			  --snp_nbhd $(FACETS_SNP_NBHD) \
			  $(if $(facets_pre_cval.$1),--pre_cval $(facets_pre_cval.$1),--pre_cval $(FACETS_PRE_CVAL)) \
			  --ndepth_max $(FACETS_NDEPTH_MAX) \
			  --use_emcncf2 \
			  $(if $(facets_diplogr.$1),--diplogr $(facets_diplogr.$1)) \
			  $(if $(facets_purity.$1),--purity $(facets_purity.$1))
SNP_PILEUP = snp-pileup
SNP_PILEUP_OPTS = -A --min-map-quality=15 --min-base-quality=15 --gzip --max-depth=150000

fastaccess : $(foreach sample,$(TUMOR_SAMPLES),fastaccess/pileup/$(sample)_A.gz fastaccess/pileup/$(sample)_B.gz)

fastaccess/vcf/targets_dbsnp_pool_A.vcf : $(POOL_A_BED)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@
	
fastaccess/vcf/targets_dbsnp_pool_B.vcf : $(POOL_B_BED)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@

define snp-pileup-tumor-normal
fastaccess/pileup/$1_A.gz : bam/$1.bam bam/$2.bam facets/vcf/targets_dbsnp_pool_A.vcf
	$$(call RUN,-c -s 8G -m 20G,"rm -f $$@ && $$(SNP_PILEUP) $$(SNP_PILEUP_OPTS) $$(<<<) $$@ $$(<<) $$(<)")
	
fastaccess/pileup/$1_B.gz : bam/$1.bam bam/$2.bam facets/vcf/targets_dbsnp_pool_B.vcf
	$$(call RUN,-c -s 8G -m 20G,"rm -f $$@ && $$(SNP_PILEUP) $$(SNP_PILEUP_OPTS) $$(<<<) $$@ $$(<<) $$(<)")

endef
 $(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


include modules/variant_callers/gatk.mk
include modules/bam_tools/processBam.mk

.PHONY : $(PHONY)
