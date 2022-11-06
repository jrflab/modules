include modules/Makefile.inc

LOGDIR ?= log/pyclone_vi.$(NOW)

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 50000'

pyclone : $(foreach sample,$(TUMOR_SAMPLES),pyclone_vi/$(sample)/$(sample).vcf) \
	  $(foreach sample,$(TUMOR_SAMPLES),pyclone_vi/$(sample)/$(sample).txt) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_vi/$(set)/$(set).tsv) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_vi/$(set)/$(set).hd5) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_vi/$(set)/summary/by_loci.txt) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_vi/$(set)/summary/scatter_by_sample.pdf) \
	  $(foreach set,$(SAMPLE_SETS),pyclone_vi/$(set)/summary/heatmap_by_sample.pdf)


define r-sufam
pyclone_vi/$1/$1.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					 --option 1 \
					 --sample_set '$(set.$1)' \
					 --normal_sample '$(normal.$1)' \
					 --input_file $$(<) \
					 --output_file $$(@)")
					 
pyclone_vi/$1/$1.txt : pyclone_vi/$1/$1.vcf bam/$1.bam
	$$(call RUN,-c -n 1 -s 2G -m 3G -v $(SUFAM_ENV),"set -o pipefail && \
					 		 sufam \
							 --sample_name $1 \
							 $$(SUFAM_OPTS) \
							 $$(REF_FASTA) \
							 $$(<) \
							 $$(<<) \
							 > $$(@)")
							 
pyclone_vi/$1/$1.maf : pyclone_vi/$1/$1.vcf
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
		
define r-pyclone
pyclone_vi/$1/$1.tsv : $(foreach sample,$(TUMOR_SAMPLES),pyclone_vi/$(sample)/$(sample).txt)
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_vi.R \
							   --option 1 \
							   --sample_set $1 \
							   --normal_sample '$(normal.$1)' \
							   --output_file $$(@)")
							   
pyclone_vi/$1/$1.hd5 : pyclone_vi/$1/$1.tsv
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(PYCLONE_ENV) -w 72:00:00,"set -o pipefail && \
							   		 pyclone-vi fit \
									 --in-file $$(<) \
									 --out-file $$(@) \
									 --num-clusters 10 \
									 --density beta-binomial \
									 --num-grid-points 100 \
									 --max-iters 1000000 \
									 --mix-weight-prior 1 \
									 --precision 500 \
									 --num-restarts 100")
									 
pyclone_vi/$1/summary/by_loci.txt : pyclone_vi/$1/$1.hd5
	$$(call RUN,-c -n 1 -s 8G -m 12G -v $(PYCLONE_ENV),"set -o pipefail && \
							   pyclone-vi write-results-file \
							   --in-file $$(<) \
							   --out-file $$(@)")
							     
pyclone_vi/$1/summary/scatter_by_sample.pdf : pyclone_vi/$1/summary/by_loci.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_vi.R \
							   --option 2 \
							   --sample_set '$(tumors.$1)' \
							   --input_file $$(<) \
							   --output_file $$(@)")
							   
pyclone_vi/$1/summary/heatmap_by_sample.pdf : pyclone_vi/$1/summary/by_loci.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G -v $(PYCLONE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/pyclone_vi.R \
							   --option 3 \
							   --sample_set '$(tumors.$1)' \
							   --input_file $$(<) \
							   --output_file $$(@)")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call r-pyclone,$(set))))
		
..DUMMY := $(shell mkdir -p version; \
	     R --version > version/pyclone_vi.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pyclone
