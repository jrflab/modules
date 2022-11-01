include modules/Makefile.inc

LOGDIR ?= log/pyclone_multi_sample.$(NOW)

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000 --ff UNMAP,SECONDARY,QCFAIL'

pyclone : $(foreach sample,$(TUMOR_SAMPLES),pyclone/$(sample).vcf) \
	  $(foreach sample,$(TUMOR_SAMPLES),pyclone/$(sample).txt) \
	  $(foreach sample,$(TUMOR_SAMPLES),pyclone/$(sample).maf)


define sufam-gt
sufam/$1.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					 --option 1 \
					 --sample_set '$(set.$1)' \
					 --normal_sample '$(normal.$1)' \
					 --input_file $$(<) \
					 --output_file $$(@)")
					 
sufam/$1.txt : sufam/$1.vcf bam/$1.bam
	$$(call RUN,-c -n 1 -s 2G -m 3G -v $(SUFAM_ENV),"set -o pipefail && \
					 		 sufam \
							 --sample_name $1 \
							 $$(SUFAM_OPTS) \
							 $$(REF_FASTA) \
							 $$(<) \
							 $$(<<) \
							 > $$(@)")

sufam/$1.maf : sufam/$1.vcf
	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(VEP_ENV),"set -o pipefail && \
							$$(VCF2MAF) \
							--input-vcf $$< \
							--tumor-id $1 \
							--filter-vcf $$(EXAC_NONTCGA) \
							--ref-fasta $$(REF_FASTA) \
							--vep-path $$(VEP_PATH) \
							--vep-data $$(VEP_DATA) \
							--tmp-dir `mktemp -d` \
							--output-maf $$(@)")


endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call sufam-gt,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
	     R --version > version/pyclone_multi_sample.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pyclone