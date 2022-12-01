include modules/Makefile.inc

LOGDIR ?= log/jags_mcmc.$(NOW)

noise_model : $(foreach sample,$(SAMPLES), \
			$(foreach barcode,$(BARCODES),hbm/$(sample)/mcmc/$(barcode).RData)) \
	      $(foreach sample,$(SAMPLES), \
			$(foreach barcode,$(BARCODES),hbm/$(sample)/params/$(barcode).RData)) \
	      $(foreach sample,$(SAMPLES),hbm/$(sample)/__LambdaP_dP__.RData)
	       

define jags-mcmc
hbm/$1/mcmc/$2.RData : summary/$1/sum_alt/$2.txt.gz vcf/$1.txt vcf/MSKCC_Weigelt_Mission_Bio_11132018.txt
	$$(call RUN,-n 1 -s 8G -m 16G -v $(JAGS_ENV) -w 72:00:00,"set -o pipefail && \
								  $(RSCRIPT) $(SCRIPTS_DIR)/mission_bio/hierarchical_bayes.R \
								  --option 1 \
								  --snp_file $$(<<) \
								  --context_file $$(<<<) \
								  --sample_name $1 \
								  --bar_code $2")

hbm/$1/params/$2.RData : hbm/$1/mcmc/$2.RData
	$$(call RUN,-n 1 -s 8G -m 16G -v $(JAGS_ENV) -w 72:00:00,"set -o pipefail && \
								  $(RSCRIPT) $(SCRIPTS_DIR)/mission_bio/hierarchical_bayes.R \
								  --option 2 \
								  --sample_name $1 \
								  --bar_code $2")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach barcode,$(BARCODES), \
		$(eval $(call jags-mcmc,$(sample),$(barcode)))))
		
define summarize-mcmc
hbm/$1/__LambdaP_dP__.RData : $(foreach barcode,$(BARCODES),hbm/$1/params/$(barcode).RData)
	$$(call RUN,-n 1 -s 24G -m 48G -v $(JAGS_ENV),"set -o pipefail && \
						       $(RSCRIPT) $(SCRIPTS_DIR)/mission_bio/hierarchical_bayes.R \
						       --option 3 \
						       --sample_name $1 \
						       --bar_code '$(BARCODES)'")

endef
$(foreach sample,$(SAMPLES), \
	$(eval $(call summarize-mcmc,$(sample))))


..DUMMY := $(shell mkdir -p version; \
	     $(JAGS_ENV)/bin/R --version > version/jags_mcmc.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: noise_model
