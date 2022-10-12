include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/annotate_maf_vcf.$(NOW)

annotate_maf_vcf : $(foreach sample,$(SAMPLES),vcf/$(sample).vcf)

define annotate-maf-vcf
vcf/%.vcf : vcf/%.maf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(VCF2MAF_ENV),"set -o pipefail && \
							   $$(MAF2VCF) \
							   --input-maf $$(<) \
							   --output-dir vcf \
							   --output-vcf $$(@) \
							   --ref-fasta $$(HOME)/share/lib/resource_files/VEP/GRCh37/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call annotate-maf-vcf,$(sample))))

							  
..DUMMY := $(shell mkdir -p version; \
	     source $(VCF2MAF_ENV)/bin/activate $(VCF2MAF_ENV) && $(MAF2VCF) --man >> version/annotate_maf_vcf.txt)
.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: annotate_maf_vcf
