include modules/Makefile.inc

LOGDIR ?= log/sufam_gt.$(NOW)

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000 --ff UNMAP,SECONDARY,QCFAIL'

sufam_gt : $(foreach sample,$(TUMOR_SAMPLES),sufam/$(sample).vcf) \
	   $(foreach sample,$(TUMOR_SAMPLES),sufam/$(sample).txt) \
	   $(foreach sample,$(TUMOR_SAMPLES),sufam/$(sample).maf)
#	   $(foreach set,$(SAMPLE_SETS),sufam/$(set).maf) \
#	   $(foreach set,$(SAMPLE_SETS),sufam/$(set)_ann.maf) \
#	   sufam/mutation_summary.maf \
#	   sufam/mutation_summary_ft.maf

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
		
define combine-maf
sufam/$1.maf : $(foreach sample,$(TUMOR_SAMPLES),sufam/$(sample).txt) $(foreach sample,$(TUMOR_SAMPLES),sufam/$(sample).maf)
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					 --option 2 \
					 --sample_set '$(set.$1)' \
					 --normal_sample '$(normal.$1)' \
					 --output_file $$(@)")
					 
sufam/$1_ft.maf : sufam/$1.maf
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					 --option 3 \
					 --input_file $$(<) \
					 --output_file $$(@)")

					 
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-maf,$(set))))
		

sufam/mutation_summary.maf : summary/tsv/all.tsv $(foreach set,$(SAMPLE_SETS),sufam/$(set).maf)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					  --option 4 \
					  --sample_set '$(SAMPLE_SETS)' \
					  --input_file $(<) \
					  --output_file $(@)")


sufam/mutation_summary_ann.maf : summary/tsv/all.tsv $(foreach set,$(SAMPLE_SETS),sufam/$(set)_ft.maf)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					  --option 5 \
					  --sample_set '$(SAMPLE_SETS)' \
					  --input_file $(<) \
					  --output_file $(@)")

..DUMMY := $(shell mkdir -p version; \
	     R --version > version/sufam_gt.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: sufam_gt
