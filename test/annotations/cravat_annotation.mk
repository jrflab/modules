include modules/Makefile.inc

LOGDIR ?= log/cravat.$(NOW)
PHONY += cravat

cravat : $(foreach sample,$(SAMPLES),cravat/$(sample).vcf)

DEFAULT_ENV = $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.6
CRAVAT_ENV = $(HOME)/share/usr/anaconda-envs/open-cravat

define cravat-annotation
cravat/%.vcf : vcf_ann/%.gatk_snps.vcf vcf_ann/%.gatk_indels.vcf
	$$(call RUN,-c -s 6G -m 8G -v $$(DEFAULT_ENV),"$(RSCRIPT) modules/test/annotation/combine_vcf.R --sample_name $$(<)")

#cravat/%.xlsx : cravat/%.vcf
#	$$(call RUN,-c -s 6G -m 8G -v $$(DEFAULT_ENV),"source activate $$(CRAVAT_ENV) && \
																 ")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cravat-annotation,$(sample))))
		
.PHONY: $(PHONY)

