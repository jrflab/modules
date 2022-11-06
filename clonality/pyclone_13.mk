include modules/Makefile.inc

LOGDIR ?= log/pyclone_13.$(NOW)

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 50000'

MCMC_ITER = 10000
MCMC_BURNIN = 2000
MCMC_THIN = 1

pyclone : $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).vcf) \
	  $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).txt) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/taskcomplete) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/config.yaml) \
	  $(foreach set,$(SAMPLE_SETS), \
	  		$(foreach sample,$(tumors.$(set)),pyclone_13/$(set)/$(sample).yaml)) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/trace/alpha.tsv.bz2) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/summary/by_clusters.txt) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/summary/by_loci.txt) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/summary/scatter_by_sample.pdf) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/summary/heatmap_by_sample.pdf)


define r-sufam
pyclone_13/$1/$1.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					 --option 1 \
					 --sample_set '$(set.$1)' \
					 --normal_sample '$(normal.$1)' \
					 --input_file $$(<) \
					 --output_file $$(@)")
					 
pyclone_13/$1/$1.txt : pyclone_13/$1/$1.vcf bam/$1.bam
	$$(call RUN,-c -n 1 -s 2G -m 3G -v $(SUFAM_ENV),"set -o pipefail && \
					 		 sufam \
							 --sample_name $1 \
							 $$(SUFAM_OPTS) \
							 $$(REF_FASTA) \
							 $$(<) \
							 $$(<<) \
							 > $$(@)")
							 
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call r-sufam,$(sample))))
		
define r-pyclone-input
pyclone_13/$1/taskcomplete : $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).txt)
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_13.R \
							   --option 1 \
							   --sample_set $1 \
							   --normal_sample '$(normal.$1)' && \
							   echo 'success' > $$(@)")
							   
pyclone_13/$1/config.yaml : $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).txt)
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_13.R \
							   --option 2 \
							   --sample_set $1 \
							   --normal_sample '$(normal.$1)' \
							   --output_file $$(@) \
							   --num_iter $$(MCMC_ITER)")
							   
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call r-pyclone-input,$(set))))
		
define r-pyclone-build-mutations
pyclone_13/$1/$2.yaml : pyclone_13/$1/taskcomplete pyclone_13/$1/config.yaml
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(PYCLONE_13_ENV),"set -o pipefail && \
							      PyClone build_mutations_file \
							      --in_file pyclone_13/$1/$2.tsv \
							      --out_file $$(@) \
							      --prior total_copy_number")
							   
endef
$(foreach set,$(SAMPLE_SETS),\
	$(foreach sample,$(tumors.$(set)),\
		$(eval $(call r-pyclone-build-mutations,$(set),$(sample)))))
		
define r-pyclone-run-analysis
pyclone_13/$1/trace/alpha.tsv.bz2 : $(foreach sample,$(tumors.$1),pyclone_13/$1/$(sample).yaml) pyclone_13/$1/config.yaml
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(PYCLONE_13_ENV) -w 72:00:00,"set -o pipefail && \
									   PyClone run_analysis \
									   --config_file pyclone_13/$1/config.yaml")
									   
pyclone_13/$1/summary/by_clusters.txt : pyclone_13/$1/trace/alpha.tsv.bz2 pyclone_13/$1/config.yaml
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(PYCLONE_13_ENV),"set -o pipefail && \
							       PyClone build_table \
							       --config_file $$(<<) \
							       --out_file $$(@) \
							       --table_type cluster \
							       --burnin $$(MCMC_BURNIN) \
							       --thin $$(MCMC_THIN)")
							       
pyclone_13/$1/summary/by_loci.txt : pyclone_13/$1/trace/alpha.tsv.bz2 pyclone_13/$1/config.yaml
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(PYCLONE_13_ENV),"set -o pipefail && \
							       PyClone build_table \
							       --config_file $$(<<) \
							       --out_file $$(@) \
							       --table_type loci \
							       --burnin $$(MCMC_BURNIN) \
							       --thin $$(MCMC_THIN)")
							       
pyclone_13/$1/summary/scatter_by_sample.pdf : pyclone_13/$1/summary/by_loci.txt pyclone_13/$1/summary/by_clusters.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_13.R \
							   --option 3 \
							   --sample_set '$(tumors.$1)' \
							   --input_file $$(<) \
							   --output_file $$(@)")
							   
pyclone_13/$1/summary/heatmap_by_sample.pdf : pyclone_13/$1/summary/by_loci.txt pyclone_13/$1/summary/by_clusters.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_13.R \
							   --option 4 \
							   --sample_set '$(tumors.$1)' \
							   --input_file $$(<) \
							   --output_file $$(@)")

							   
endef
$(foreach set,$(SAMPLE_SETS),\
	$(eval $(call r-pyclone-run-analysis,$(set))))
		

..DUMMY := $(shell mkdir -p version; \
	     R --version > version/pyclone_13.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pyclone
