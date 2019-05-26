include modules/Makefile.inc

LOGDIR ?= log/fastaccess.$(NOW)
PHONY += fastaccess fastaccess/vcf fastaccess/snp fastaccess/cnr fastaccess/plots

FACETS_ENV = $(HOME)/share/usr/anaconda-envs/facets-0.5.6/
RUN_FACETS = $(RSCRIPT) modules/test/copy_number/fastaccess.R
POOL_A_BED = ~/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
POOL_B_BED = ~/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
SNP_PILEUP = snp-pileup
SNP_PILEUP_OPTS = -A --min-map-quality=15 --min-base-quality=15 --gzip --max-depth=150000

fastaccess : $(foreach sample,$(TUMOR_SAMPLES),fastaccess/snp/$(sample)_A.gz fastaccess/snp/$(sample)_B.gz fastaccess/cnr/$(sample).txt)

fastaccess/vcf/targets_dbsnp_pool_A.vcf : $(POOL_A_BED)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@
	
fastaccess/vcf/targets_dbsnp_pool_B.vcf : $(POOL_B_BED)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@

define snp-pileup-tumor-normal
fastaccess/snp/$1_A.gz : bam/$1.bam bam/$2.bam fastaccess/vcf/targets_dbsnp_pool_A.vcf
	$$(call RUN,-c -s 8G -m 20G,"rm -f $$@ && $$(SNP_PILEUP) $$(SNP_PILEUP_OPTS) $$(<<<) $$@ $$(<<) $$(<)")
	
fastaccess/snp/$1_B.gz : bam/$1.bam bam/$2.bam fastaccess/vcf/targets_dbsnp_pool_B.vcf
	$$(call RUN,-c -s 8G -m 20G,"rm -f $$@ && $$(SNP_PILEUP) $$(SNP_PILEUP_OPTS) $$(<<<) $$@ $$(<<) $$(<)")

endef
 $(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

 
define fast-access-tumor
fastaccess/cnr/%.txt : fastaccess/snp/%_A.gz fastaccess/snp/%_B.gz
	$$(call RUN,-c -v $(FACETS_ENV) -s 8G -m 60G,"$(RUN_FACETS) --option 1 --pool_A $$(<) --pool_B $$(<<)")

endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call fast-access-tumor,$(sample)))) 
 

include modules/variant_callers/gatk.mk
include modules/bam_tools/processBam.mk

.PHONY : $(PHONY)
