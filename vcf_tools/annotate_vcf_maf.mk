include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/annotate_vcf_maf.$(NOW)

annotate_vcf_maf : $(foreach sample,$(SAMPLES),vcf/$(sample).maf)

define annotate-vcf-maf
vcf/%.maf : vcf/%.vcf
	$$(call RUN,-c -n 12 -s 2G -m 3G -v $(VCF2MAF_ENV) -w 72:00:00,"set -o pipefail && \
									$$(VCF2MAF) \
									--input-vcf $$(<) \
									--output-maf $$(@) \
									--tmp-dir $$(TMPDIR) \
									--tumor-id $$(*) \
									--normal-id NA \
									--vep-path $$(VCF2MAF_ENV)/bin \
									--vep-data $$(HOME)/share/lib/resource_files/VEP/GRCh37/ \
									--vep-forks 12 \
									--ref-fasta $$(HOME)/share/lib/resource_files/VEP/GRCh37/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
									--filter-vcf $$(HOME)/share/lib/resource_files/VEP/GRCh37/homo_sapiens/99_GRCh37/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
									--species homo_sapiens \
									--ncbi-build GRCh37 \
									--maf-center MSKCC && \
									$$(RM) $$(TMPDIR)/$$(*).vep.vcf")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call annotate-vcf-maf,$(sample))))

							  
..DUMMY := $(shell mkdir -p version; \
	     source $(VCF2MAF_ENV)/bin/activate $(VCF2MAF_ENV) && $(VCF2MAF) --man >> version/annotate_vcf_maf.txt)
.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: annotate_vcf_maf
