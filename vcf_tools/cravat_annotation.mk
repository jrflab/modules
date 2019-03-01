include modules/Makefile.inc

LOGDIR ?= log/cravat_annotate.$(NOW)
PHONY += cravat

cravat_annotate : $(foreach sample,$(SAMPLES),cravat/$(sample).vcf cravat/$(sample).maf cravat/$(sample).cravat.vcf cravat/$(sample).tsv cravat/$(sample).txt)

DEFAULT_ENV = $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.6
CRAVAT_ENV = $(HOME)/share/usr/anaconda-envs/open-cravat

define cravat-annotation
cravat/%.vcf : vcf_ann/%.gatk_snps.vcf vcf_ann/%.gatk_indels.vcf
	$$(call RUN,-c -s 9G -m 12G -w 7200,"$(RSCRIPT) modules/vcf_tools/combine_vcf.R --sample_name $$(*)")
	
cravat/%.maf : cravat/%.vcf
	$$(call RUN,-s 9G -m 12G -v $$(VEP_ENV) -w 7200,"$$(VCF2MAF) --input-vcf $$< --tumor-id $$(*) $$(if $$(EXAC_NONTCGA),--filter-vcf $$(EXAC_NONTCGA)) --ref-fasta $$(REF_FASTA) --vep-path $$(VEP_PATH) --vep-data $$(VEP_DATA) --tmp-dir `mktemp -d` --output-maf $$@")

cravat/%.cravat.vcf : cravat/%.vcf cravat/%.maf
	$$(call RUN,-c -s 9G -m 12G -w 7200,"$(RSCRIPT) modules/vcf_tools/filter_vcf.R --sample_name $$(*)")

cravat/%.tsv: cravat/%.cravat.vcf
	$$(call RUN,-c -s 9G -m 12G -v $$(DEFAULT_ENV) -w 7200,"source activate $$(CRAVAT_ENV) && \
												    		cravat cravat/$$(*).cravat.vcf -n $$(*) -d cravat -a clinvar cosmic dbsnp gnomad hgvs -v -l hg19 -t text")
												    
cravat/%.txt : cravat/%.tsv
	$$(call RUN,-c -s 9G -m 12G -w 7200,"$(RSCRIPT) modules/vcf_tools/summary_vcf.R --sample_name $$(*)")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cravat-annotation,$(sample))))
		
.PHONY: $(PHONY)

