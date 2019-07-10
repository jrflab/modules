include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_reference.$(NOW)
PHONY += cnvkit cnvaccess/reference

cnvaccess_reference: cnvaccess/reference/on_target_pool_A.cnr cnvaccess/reference/on_target_pool_B.cnr cnvaccess/reference/off_target_no_pool.cnr

cnvaccess/reference/on_target_pool_A.cnr : $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-A.targetcoverage.cnn) $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-A.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvaccess/cnn/normal/*.pool-A.*.cnn -f $(REF_FASTA) --no-edge -o cnvaccess/reference/on_target_pool_A.cnr")
	
cnvaccess/reference/on_target_pool_B.cnr : $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-B.targetcoverage.cnn) $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-B.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvaccess/cnn/normal/*.pool-B.*.cnn -f $(REF_FASTA) --no-edge -o cnvaccess/reference/on_target_pool_B.cnr")
	
cnvaccess/reference/off_target_no_pool.cnr : $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).no-pool.targetcoverage.cnn) $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).no-pool.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvaccess/cnn/normal/*.no-pool.*.cnn -f $(REF_FASTA) --no-edge -o cnvaccess/reference/off_target_no_pool.cnr")

		
.PHONY: $(PHONY)
