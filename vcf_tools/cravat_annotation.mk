include modules/Makefile.inc

LOGDIR ?= log/cravat_annotate.$(NOW)

cravat_annotate : $(foreach sample,$(SAMPLES),cravat/$(sample).vcf) \
		  $(foreach sample,$(SAMPLES),cravat/$(sample).maf) \
		  $(foreach sample,$(SAMPLES),cravat/$(sample).cravat.vcf) \
		  $(foreach sample,$(SAMPLES),cravat/$(sample).tsv) \
		  $(foreach sample,$(SAMPLES),cravat/$(sample).txt)

define cravat-annotation
cravat/$1.vcf : vcf_ann/$1.gatk_snps.vcf vcf_ann/$1.gatk_indels.vcf
	$$(call RUN,-c -s 9G -m 12G -w 24:00:00,"set -o pipefail && \
						 $(RSCRIPT) modules/vcf_tools/combine_vcf.R \
						 --sample_name $$(*)")
	
cravat/$1.maf : cravat/$1.vcf
	$$(call RUN,-c -s 9G -m 12G -v $(VEP_ENV) -w 24:00:00,"set -o pipefail && \
							       $$(VCF2MAF) \
							       --input-vcf $$(<) \
							       --tumor-id $1 \
							       $$(if $$(EXAC_NONTCGA),--filter-vcf $$(EXAC_NONTCGA)) \
							       --ref-fasta $$(REF_FASTA) \
							       --vep-path $$(VEP_PATH) \
							       --vep-data $$(VEP_DATA) \
							       --tmp-dir `mktemp -d` \
							       --output-maf $$(@)")

cravat/$1.cravat.vcf : cravat/$1.vcf cravat/$1.maf
	$$(call RUN,-c -s 9G -m 12G -w 24:00:00,"set -o pipefail && \
						 $(RSCRIPT) modules/vcf_tools/filter_vcf.R \
						 --sample_name $1")

cravat/$1.tsv: cravat/$1.cravat.vcf
	$$(call RUN,-c -s 9G -m 12G -v $(CRAVAT_ENV) -w 24:00:00,"set -o pipefail && \
								  cravat $$(<) \
								  -n $1 \
								  -d cravat \
								  -a clinvar cosmic dbsnp gnomad hgvs \
								  -v \
								  -l hg19 \
								  -t text")
												    
cravat/$1.txt : cravat/$1.tsv
	$$(call RUN,-c -s 9G -m 12G -w 24:00:00,"set -o pipefail && \
						 $(RSCRIPT) modules/vcf_tools/summary_vcf.R \
						 --sample_name $1")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cravat-annotation,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
             echo "cravat" > version/cravat_annotate.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : cravat_annotate
