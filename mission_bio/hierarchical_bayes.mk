include modules/Makefile.inc

LOGDIR ?= log/hierarchical_bayes.$(NOW)

noise_model : $(foreach sample,$(SAMPLES), \
			$(foreach barcode,$(BARCODES),hbm/$(sample)/mcmc/$(barcode).RData))
	       

define hierarchical-bayes
hbm/$1/mcmc/$2.RData : summary/$1/sum_alt/$2.txt.gz vcf/$1.txt vcf/MSKCC_Weigelt_Mission_Bio_11132018.txt
	$$(call RUN,-n 1 -s 8G -m 16G -v $(JAGS_ENV) -w 24:00:00,"set -o pipefail && \
								   $(RSCRIPT) $(SCRIPTS_DIR)/mission_bio/hierarchical_bayes.R \
								   --option 1 \
								   --snp_file $$(<<) \
								   --context_file $$(<<<) \
								   --sample_name $1 \
								   --bar_code $2")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach barcode,$(BARCODES), \
		$(eval $(call hierarchical-bayes,$(sample),$(barcode)))))


..DUMMY := $(shell mkdir -p version; \
	     $(JAGS_ENV)/bin/R --version > version/hierarchical_bayes.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: noise_model
