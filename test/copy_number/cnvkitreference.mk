include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_reference.$(NOW)
PHONY += cnvkit cnvkit/reference

cnvkit_reference: cnvkit/reference/on_target_pool_A.cnr cnvkit/reference/on_target_pool_B.cnr cnvkit/reference/off_target_no_pool.cnr

cnvkit/reference/on_target_pool_A.cnr : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).A.targetcoverage.cnn) $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).A.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvkit/cnn/normal/*.A.*.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/on_target_pool_A.cnr")
	
cnvkit/reference/on_target_pool_B.cnr : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).B.targetcoverage.cnn) $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).B.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvkit/cnn/normal/*.B.*.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/on_target_pool_B.cnr")
	
cnvkit/reference/off_target_no_pool.cnr : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).C.targetcoverage.cnn) $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).C.antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvkit/cnn/normal/*.C.*.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/off_target_no_pool.cnr")

		
.PHONY: $(PHONY)
