include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/annotate_vcf_context.$(NOW)

annotate_vcf_context : $(foreach sample,$(SAMPLES),vcf/$(sample).txt)

define annotate-vcf-context
vcf/%.txt : vcf/%.vcf
	$$(call RUN,-c -s 4G -m 8G -v $(ANNOTATION_ENV),"set -o pipefail && \
							 $$(RSCRIPT) $$(SCRIPTS_DIR)/vcf_tools/annotate_vcf_context.R \
                                                   	 --file_in $$(<) \
                                                   	 --file_out $$(@) \
                                                   	 --ensembl_gene $$(ENSEMBL)")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call annotate-vcf-context,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     R --version >> version/annotate_vcf_context.txt)
.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: annotate_vcf_context
