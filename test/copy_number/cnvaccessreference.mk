include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_reference.$(NOW)
PHONY += cnvkit cnvaccess/reference

ACCESS_REF_FILE_A ?= ~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-A.cnr
ACCESS_REF_FILE_B ?= ~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-B.cnr
ACCESS_REF_FILE_OFF ?= ~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-noprobe.cnr

USE_REF ?= true

cnvaccess_reference: cnvaccess/reference/on_target_pool_A.cnr cnvaccess/reference/on_target_pool_B.cnr cnvaccess/reference/off_target_no_pool.cnr

ifeq ($(USE_REF),true)

cnvaccess/reference/on_target_pool_A.cnr : $(ACCESS_REF_FILE_A)
	$(call RUN,-n 1 -s 1G -m 2G,"cp $$(ACCESS_REF_FILE_A) cnvaccess/reference/on_target_pool_A.cnr")
	
cnvaccess/reference/on_target_pool_B.cnr : $(ACCESS_REF_FILE_B)
	$(call RUN,-n 1 -s 1G -m 2G,"cp $$(ACCESS_REF_FILE_B) cnvaccess/reference/on_target_pool_B.cnr")
	
cnvaccess/reference/off_target_no_pool.cnr : $(ACCESS_REF_FILE_OFF)
	$(call RUN,-n 1 -s 1G -m 2G,"cp $$(ACCESS_REF_FILE_OFF) cnvaccess/reference/off_target_no_pool.cnr")

else

cnvaccess/reference/on_target_pool_A.cnr : $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-A.targetcoverage.cnn) $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-A.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvaccess/cnn/normal/*.pool-A.*.cnn -f $(REF_FASTA) --no-edge -o cnvaccess/reference/on_target_pool_A.cnr")
	
cnvaccess/reference/on_target_pool_B.cnr : $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-B.targetcoverage.cnn) $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).pool-B.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvaccess/cnn/normal/*.pool-B.*.cnn -f $(REF_FASTA) --no-edge -o cnvaccess/reference/on_target_pool_B.cnr")
	
cnvaccess/reference/off_target_no_pool.cnr : $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).no-pool.targetcoverage.cnn) $(wildcard cnvaccess/cnn/normal/$(NORMAL_SAMPLES).no-pool.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvaccess/cnn/normal/*.no-pool.*.cnn -f $(REF_FASTA) --no-edge -o cnvaccess/reference/off_target_no_pool.cnr")

endif
		
.PHONY: $(PHONY)
