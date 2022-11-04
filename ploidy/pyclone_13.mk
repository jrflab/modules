include modules/Makefile.inc

LOGDIR ?= log/pyclone_13.$(NOW)

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 50000'

pyclone : $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).vcf) \
	  $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).txt) \
	  $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).maf) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/$(set).taskcomplete) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/config.yaml) \
	  $(foreach set,$(SAMPLE_SETS), \
	  		$(foreach sample,$(tumors.$(set)),pyclone_13/$(set)/$(sample).yaml)) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_13/$(set)/trace/alpha.tsv.bz2)
#	  $(foreach set,$(SAMPLE_SETS),pyclone_vi/$(set)/$(set)__PS__.pdf) \
#	  $(foreach set,$(SAMPLE_SETS),pyclone_vi/$(set)/$(set)__HM__.pdf)


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
							 
pyclone_13/$1/$1.maf : pyclone_13/$1/$1.vcf
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
		$(eval $(call r-sufam,$(sample))))
		
define r-pyclone-input
pyclone_13/$1/$1.taskcomplete : $(foreach sample,$(TUMOR_SAMPLES),pyclone_13/$(sample)/$(sample).txt)
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_13.R \
							   --option 1 \
							   --sample_set $1 \
							   --normal_sample '$(normal.$1)' && \
							   echo 'taskcomplete' > $$(@)")
							   
pyclone_13/$1/config.yaml : pyclone_13/$1/$1.taskcomplete
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_13.R \
							   --option 2 \
							   --sample_set $1 \
							   --normal_sample '$(normal.$1)' \
							   --output_file $$(@)")
							   
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call r-pyclone-input,$(set))))
		
define r-pyclone-build-mutations
pyclone_13/$1/$2.yaml : pyclone_13/$1/$1.taskcomplete pyclone_13/$1/config.yaml
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
pyclone_13/$1/trace/alpha.tsv.bz2 : $(foreach sample,$(tumors.$1),pyclone_13/$1/$(sample).yaml)
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(PYCLONE_13_ENV) -w 72:00:00,"set -o pipefail && \
									   PyClone run_analysis \
									   --config_file pyclone_13/$1/config.yaml")
							   
endef
$(foreach set,$(SAMPLE_SETS),\
	$(eval $(call r-pyclone-run-analysis,$(set))))
		

..DUMMY := $(shell mkdir -p version; \
	     R --version > version/pyclone_13.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pyclone
