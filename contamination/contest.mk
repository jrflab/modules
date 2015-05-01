# This module runs ContEst on snp vcf files from gatk
# Author: inodb

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/contest.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: contest

contest : contest/snp_vcf/all_contamination.txt $(foreach sample,$(SAMPLES),contest/snp_vcf/$(sample)_contamination.txt)

CONTEST_JAR ?= ~/share/usr/contest-1.0.24530-bin/ContEst.jar
CONTEST_REF_FA ?= ~/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
CONTEST_REF_VCF ?= ~/share/reference/hg19_population_stratified_af_hapmap_3.3.vcf
CONTEST_MEM = $(JAVA) -Xmx$(1) -jar $(CONTEST_JAR) -R $(CONTEST_REF_FA) -B:pop$(,)vcf $(CONTEST_REF_VCF) -T Contamination -BTI genotypes

# ContEst on gatk snp_vcf folder
contest/snp_vcf/%_contamination.txt : bam/%.bam snp_vcf/%.snps.vcf
	$(call LSCRIPT_MEM,7G,8G,"$(MKDIR) contest/snp_vcf; $(call CONTEST_MEM,2G) -I $(<) -B:genotypes$(,)vcf $(<<) -o $(@)")

contest/snp_vcf/all_contamination.txt : $(foreach sample,$(SAMPLES),contest/snp_vcf/$(sample)_contamination.txt)
	( \
		head -1 $< | paste <(echo sample) -; \
		for s in $(^); do \
			grep "META	" $$s | paste <(echo `basename $$s _contamination.txt`) -; \
		done | sort -rnk5,5; \
	) > $@
